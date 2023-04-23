//! Prototype streamlined storage of haploid genomes
//! where we also count the number of occurrences
//! of each genome.
//!
//! The goal is to prevent unnecessary copies of
//! parent genomes to offspring.
//!
//! This started as a copy/paste of genome_array.rs

use rand::prelude::Rng;
use rand::SeedableRng;

use forrustts::genetics::{GenerateBreakpoints, GeneticMap};
use forrustts::prelude::*;

use crate::common::generate_mutations;
use crate::common::generate_offspring_genome;
use crate::common::set_fixation_counts_to_zero;
use crate::common::DiploidGenome;
use crate::common::Mutation;
use crate::common::MutationRange;
use crate::common::ParentalGenome;
use crate::common::SimParams;

#[derive(Debug, Default)]
struct HaploidGenomes {
    genome_spans: Vec<MutationRange>,
    // Number of occurrences of each span
    counts: Vec<u32>,
    mutations: Vec<u32>,
}

impl HaploidGenomes {
    fn get_genome(&self, genome: usize) -> ParentalGenome {
        if genome != usize::MAX {
            let index_range = self.genome_spans[genome];
            let mutations = &self.mutations[index_range.start..index_range.stop];
            ParentalGenome {
                mutations,
                current_mutation_index: 0,
                genome,
            }
        } else {
            ParentalGenome {
                mutations: &[],
                current_mutation_index: 0,
                genome,
            }
        }
    }

    fn add_range(&mut self, range: MutationRange) -> usize {
        if range.start == range.stop {
            usize::MAX
        } else {
            let rv = self.genome_spans.len();
            self.genome_spans.push(range);
            self.counts.push(1);
            rv
        }
    }

    fn clear(&mut self) {
        self.mutations.clear();
        self.genome_spans.clear();
        self.counts.clear();
    }
}

// Implementation is specific to "diploid",
// so not an associated fn of HaploidGenomes
fn get_parental_genomes(
    genomes: &HaploidGenomes,
    parent: DiploidGenome,
) -> (ParentalGenome, ParentalGenome) {
    (
        genomes.get_genome(parent.first),
        genomes.get_genome(parent.second),
    )
}

/// When will the borrow checker hate this?
pub struct DiploidPopulation {
    genomes: HaploidGenomes,
    individuals: Vec<DiploidGenome>,
    mutations: Vec<Mutation>,
    mutation_counts: Vec<u32>,
}

impl DiploidPopulation {
    pub fn new(size: u32) -> Option<Self> {
        if size > 0 {
            let genomes = HaploidGenomes::default();

            // Now, everyone starts with a single "empty"
            // genome
            let individuals = vec![DiploidGenome::new(usize::MAX, usize::MAX); size as usize];

            Some(Self {
                genomes,
                individuals,
                mutations: vec![],
                mutation_counts: vec![],
            })
        } else {
            None
        }
    }

    /// Generate indexes of extinct mutations
    #[inline(never)]
    fn mutation_recycling(&self) -> Vec<usize> {
        self.mutation_counts
            .iter()
            .enumerate()
            .filter_map(|(index, value)| if value == &0 { Some(index) } else { None })
            .collect::<Vec<usize>>()
    }

    #[inline(never)]
    pub fn count_mutations(&mut self) {
        self.mutation_counts.fill(0);
        self.mutation_counts.resize(self.mutations.len(), 0);
        for (range, count) in self
            .genomes
            .genome_spans
            .iter()
            .zip(self.genomes.counts.iter())
        {
            self.genomes.mutations[range.start..range.stop]
                .iter()
                .for_each(|m| {
                    *unsafe { self.mutation_counts.get_unchecked_mut(*m as usize) } += count
                });
        }
    }

    #[inline(never)]
    pub fn num_segregating_mutations(&self) -> u32 {
        let mut nseg = 0;
        for i in &self.mutation_counts {
            if i > &0 && (*i as usize) < 2 * self.individuals.len() {
                nseg += 1;
            }
        }
        nseg
    }

    pub fn sum_extant_genome_sizes(&self) -> usize {
        let mut sum = 0;
        for (index, value) in self.genomes.genome_spans.iter().enumerate() {
            sum += (value.stop - value.start) * (self.genomes.counts[index] as usize);
        }
        sum
    }
}

// This normally wouldn't be pub,
// but should be unit tested anyways.
#[inline(never)]
fn fixation_removal_check(mutation_counts: &[u32], twon: u32, output: &mut HaploidGenomes) -> bool {
    if mutation_counts.iter().any(|m| *m == twon) {
        let x = output.mutations.len();
        // SAFETY: see comments in genome_array.rs
        output
            .mutations
            .retain(|m| *unsafe { mutation_counts.get_unchecked(*m as usize) } < twon);
        let delta = x - output.mutations.len();
        assert_eq!(delta % output.genome_spans.len(), 0);

        let delta_per_genome = delta / output.genome_spans.len();

        // NOTE: could be SIMD later
        output
            .genome_spans
            .iter_mut()
            .enumerate()
            .for_each(|(i, h)| {
                let c = h.start;
                h.start -= i * delta_per_genome;
                if i > 0 {
                    assert!(c > h.start);
                }
                assert!(h.start < output.mutations.len());
                h.stop -= (i + 1) * delta_per_genome;
                assert!(h.stop >= h.start);
                assert!(
                    h.stop <= output.mutations.len(),
                    "{h:?} {}",
                    output.mutations.len()
                );
            });
        true
    } else {
        false
    }
}

#[derive(Copy, Clone, Debug)]
struct SimDistributions {
    u01: rand::distributions::Uniform<f64>,
    num_mutations: rand_distr::Poisson<f64>,
    position_generator: rand::distributions::Uniform<Position>,
    parent_picker: rand::distributions::Uniform<usize>,
}

impl SimDistributions {
    fn new(num_individuals: u32, mutation_rate: f64, genome_length: Position) -> Option<Self> {
        let parent_picker = rand::distributions::Uniform::<usize>::new(0, num_individuals as usize);
        let u01 = rand::distributions::Uniform::new(0., 1.);
        let num_mutations = rand_distr::Poisson::<f64>::new(mutation_rate).ok()?;
        let position_generator =
            rand::distributions::Uniform::<Position>::new(Position::new_valid(0), genome_length);
        Some(Self {
            u01,
            num_mutations,
            position_generator,
            parent_picker,
        })
    }
}

fn generate_offspring_details(
    parent: usize,
    generation: u32,
    distributions: SimDistributions,
    rng: &mut rand::rngs::StdRng,
    genetic_map: &mut GeneticMap,
    queue: &mut Vec<usize>,
    pop: &mut DiploidPopulation,
    offspring_genomes: &mut HaploidGenomes,
    output_genome_map: &mut [usize],
) -> usize {
    let mutations = generate_mutations(
        generation,
        distributions.num_mutations,
        distributions.position_generator,
        queue,
        &mut pop.mutations,
        rng,
    );

    genetic_map.generate_breakpoints(rng);
    pop.mutation_counts
        .resize(pop.mutation_counts.len() + mutations.len(), 0);

    let genomes = get_parental_genomes(&pop.genomes, pop.individuals[parent]);
    let genomes = if rng.sample(distributions.u01) < 0.5 {
        genomes
    } else {
        (genomes.1, genomes.0)
    };
    if mutations.is_empty() && genetic_map.is_empty() {
        if genomes.0.genome != usize::MAX {
            assert!(genomes.0.genome < output_genome_map.len());
            if output_genome_map[genomes.0.genome] != usize::MAX {
                offspring_genomes.counts[output_genome_map[genomes.0.genome]] += 1;
                output_genome_map[genomes.0.genome]
            } else {
                output_genome_map[genomes.0.genome] = offspring_genomes.genome_spans.len();
                let start = offspring_genomes.mutations.len();
                offspring_genomes
                    .mutations
                    .extend_from_slice(genomes.0.mutations);
                let stop = offspring_genomes.mutations.len();
                offspring_genomes
                    .genome_spans
                    .push(MutationRange { start, stop });
                offspring_genomes.counts.push(1);
                offspring_genomes.genome_spans.len() - 1
            }
        } else {
            usize::MAX
        }
    } else {
        let range = generate_offspring_genome(
            genomes,
            &pop.mutations,
            mutations,
            genetic_map.breakpoints(),
            &mut offspring_genomes.mutations,
        );
        offspring_genomes.add_range(range)
    }
}

fn generate_offspring(
    generation: u32,
    distributions: SimDistributions,
    rng: &mut rand::rngs::StdRng,
    genetic_map: &mut GeneticMap,
    queue: &mut Vec<usize>,
    pop: &mut DiploidPopulation,
    offspring_genomes: &mut HaploidGenomes,
    output_genome_map: &mut [usize],
) -> DiploidGenome {
    let parent1 = rng.sample(distributions.parent_picker);
    let parent2 = rng.sample(distributions.parent_picker);

    let first = generate_offspring_details(
        parent1,
        generation,
        distributions,
        rng,
        genetic_map,
        queue,
        pop,
        offspring_genomes,
        output_genome_map,
    );
    let second = generate_offspring_details(
        parent2,
        generation,
        distributions,
        rng,
        genetic_map,
        queue,
        pop,
        offspring_genomes,
        output_genome_map,
    );
    DiploidGenome { first, second }
}

// A proper implementation
// would be generic over "generating mutations"
#[inline(never)]
pub fn evolve_pop(params: SimParams, genetic_map: GeneticMap) -> Option<DiploidPopulation> {
    let params = params.validate()?;
    let mut pop = DiploidPopulation::new(params.num_individuals)?;

    let mut rng = rand::rngs::StdRng::seed_from_u64(params.seed);

    let distributions = SimDistributions::new(
        params.num_individuals,
        params.mutation_rate,
        Position::new_valid(1000000),
    )?;

    let mut genetic_map = genetic_map;

    //let mut parent_haplotype_map = vec![];
    let mut offspring_genomes = HaploidGenomes::default();
    offspring_genomes.mutations.reserve(1000);
    let mut offspring = vec![];
    let mut output_genome_map = vec![];
    for generation in 0..params.num_generations {
        let mut queue = pop.mutation_recycling();
        output_genome_map.fill(usize::MAX);
        output_genome_map.resize(pop.genomes.genome_spans.len(), usize::MAX);
        for _ in 0..params.num_individuals {
            let offspring_genome = generate_offspring(
                generation,
                distributions,
                &mut rng,
                &mut genetic_map,
                &mut queue,
                &mut pop,
                &mut offspring_genomes,
                &mut output_genome_map,
            );
            offspring.push(offspring_genome);
        }
        std::mem::swap(&mut pop.genomes, &mut offspring_genomes);
        offspring_genomes.clear();

        std::mem::swap(&mut pop.individuals, &mut offspring);
        offspring.clear();
        pop.count_mutations();

        if fixation_removal_check(
            &pop.mutation_counts,
            2 * params.num_individuals,
            &mut pop.genomes,
        ) {
            set_fixation_counts_to_zero(2 * params.num_individuals, &mut pop.mutation_counts);
        };
    }
    Some(pop)
}

#[cfg(test)]
mod tests_to_delete {
    #[test]
    fn test_slice_behavior() {
        let v = Vec::<i32>::new();
        let s = &v[0..];
        assert!(s.is_empty());
        let ss = &s[0..];
        assert!(ss.is_empty());
    }

    #[test]
    fn test_swap_slices() {
        let v = vec![1, 2, 3, 4];
        let mut s0 = &v[0..2];
        let mut s1 = &v[2..];
        std::mem::swap(&mut s0, &mut s1);
        assert_eq!(s0, &[3, 4]);
        assert_eq!(s1, &[1, 2]);
    }

    #[test]
    fn test_size() {
        // This is one idea of production code
        // could employ for "haplotype keys".
        struct X(std::num::NonZeroUsize);
        assert_eq!(std::mem::size_of::<X>(), std::mem::size_of::<usize>());
        assert_eq!(
            std::mem::size_of::<Option<X>>(),
            std::mem::size_of::<usize>()
        );
    }
}
