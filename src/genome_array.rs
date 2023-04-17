//! Prototype streamlined storage of haploid genomes

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
            rv
        }
    }

    fn clear(&mut self) {
        self.mutations.clear();
        self.genome_spans.clear();
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
        self.genomes
            .mutations
            .iter()
            .for_each(|m| self.mutation_counts[*m as usize] += 1);
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
        self.genomes.mutations.len()
    }
}

// fn generate_offspring_genome2(
//     parent: DiploidGenome,
//     parent_haplotypes: &Haplotypes,
//     mutations: &[Mutation],
//     new_mutations: Vec<usize>,
//     breakpoints: &[Breakpoint],
//     offspring_haplotypes: &mut Haplotypes,
//     parent_haplotype_map: &mut [usize],
//     rng: &mut rand::rngs::StdRng,
// ) -> usize {
//     let u01 = rand::distributions::Uniform::new(0., 1.);
//     let (mut first_genome, mut second_genome) = get_parental_genomes(parent_haplotypes, parent);
//     if rng.sample(u01) < 0.5 {
//         std::mem::swap(&mut first_genome, &mut second_genome);
//     }
//     let mut rv = usize::MAX;
//     if breakpoints.is_empty() && new_mutations.is_empty() {
//         // Then we have aready processed this genome
//         if parent_haplotype_map[first_genome.genome] != usize::MAX {
//             return parent_haplotype_map[first_genome.genome];
//         }
//     }
//     let start = offspring_haplotypes.mutations.len();
//     let nm = new_mutations.len();
//     for m in new_mutations.iter() {
//         let n = first_genome.mutations[first_genome.current_mutation_index..]
//             .iter()
//             .take_while(|mutation| mutations[**mutation].position() < mutations[*m].position())
//             .inspect(|x| offspring_haplotypes.mutations.push(**x))
//             .count();
//         offspring_haplotypes.mutations.push(*m);
//         first_genome.current_mutation_index += n;
//     }
//     first_genome.mutations[first_genome.current_mutation_index..]
//         .iter()
//         .for_each(|m| offspring_haplotypes.mutations.push(*m));
//     let stop = offspring_haplotypes.mutations.len();
//     if stop > start {
//         rv = offspring_haplotypes.haplotypes.len();
//         offspring_haplotypes
//             .haplotypes
//             .push(MutationRange { start, stop });
//     }
//     if nm == 0 && breakpoints.is_empty() {
//         assert_eq!(parent_haplotype_map[first_genome.genome], usize::MAX);
//         parent_haplotype_map[first_genome.genome] = rv;
//     }
//     assert_eq!(
//         stop - start,
//         nm + first_genome.mutations.len(),
//         "{:?} + {new_mutations:?} = {:?}",
//         first_genome.mutations,
//         &offspring_haplotypes.mutations[start..stop]
//     );
//     rv
// }

impl SimParams {
    pub fn validate(self) -> Option<Self> {
        if !self.mutation_rate.is_finite() || self.mutation_rate < 0.0 {
            return None;
        }
        Some(self)
    }
}

// This normally wouldn't be pub,
// but should be unit tested anyways.
#[inline(never)]
fn fixation_removal_check(mutation_counts: &[u32], twon: u32, output: &mut HaploidGenomes) -> bool {
    if mutation_counts.iter().any(|m| *m == twon) {
        let x = output.mutations.len();
        output
            .mutations
            .retain(|m| mutation_counts[*m as usize] < twon);
        let delta = x - output.mutations.len();
        assert_eq!(delta % output.genome_spans.len(), 0);

        let delta_per_genome = delta / output.genome_spans.len();
        assert_eq!(
            delta_per_genome,
            mutation_counts.iter().filter(|m| **m == twon).count()
        );

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
    let range = generate_offspring_genome(
        genomes,
        &pop.mutations,
        mutations,
        genetic_map.breakpoints(),
        &mut offspring_genomes.mutations,
    );
    if range.stop > range.start {
        offspring_genomes.add_range(range)
    } else {
        // offspring genome is empty
        usize::MAX
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
    for generation in 0..params.num_generations {
        let mut queue = pop.mutation_recycling();
        for _ in 0..params.num_individuals {
            let offspring_genome = generate_offspring(
                generation,
                distributions,
                &mut rng,
                &mut genetic_map,
                &mut queue,
                &mut pop,
                &mut offspring_genomes,
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
mod tests {
    use super::*;

    use proptest::prelude::*;

    proptest! {
        #[test]
        #[ignore]
        fn run_sim_no_recombination(seed in 0..u64::MAX) {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            let make_mutrate = rand_distr::Exp::new(1.0).unwrap();
            let mutation_rate = rng.sample(make_mutrate);
            let params = SimParams {
                seed,
                num_individuals: 100,
                num_generations: 100,
                mutation_rate,
                recrate: 0.0,
            };
            // Empty genetic map == no recombination
            let builder = forrustts::genetics::GeneticMapBuilder::default();
            let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
            let _ = evolve_pop(params, genetic_map).unwrap();
        }
    }

    proptest! {
        #[test]
        #[ignore]
        fn run_sim_with_recombination(seed in 0..u64::MAX) {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            let make_mutrate = rand_distr::Exp::new(1.0).unwrap();
            let mutation_rate = rng.sample(make_mutrate);
            let params = SimParams {
                seed,
                num_individuals: 100,
                num_generations: 100,
                mutation_rate,
                recrate: 0.0,
            };

            let genome_start = Position::new_valid(0);
            let genome_length = Position::new_valid(1000000);
            let poisson = vec![forrustts::genetics::PoissonCrossover::new(
                genome_start, genome_length, 2.0).unwrap()];
            let builder = forrustts::genetics::GeneticMapBuilder::default().extend_poisson(&poisson);

            let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
            let _ = evolve_pop(params, genetic_map).unwrap();
        }
    }
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
