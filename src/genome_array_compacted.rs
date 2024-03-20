//! Prototype streamlined storage of haploid genomes

use rand::prelude::Rng;
use rand::SeedableRng;

use forrustts::genetics::{GenerateBreakpoints, GeneticMap};
use forrustts::prelude::*;

use crate::common::generate_mutations;
use crate::common::generate_offspring_genome;
use crate::common::DiploidGenome;
use crate::common::Mutation;
use crate::common::MutationRange;
use crate::common::ParentalGenome;
use crate::common::ParentalGenomeMut;
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

    fn get_genome_mut(&mut self, genome: usize) -> ParentalGenomeMut {
        if genome != usize::MAX {
            let index_range = self.genome_spans[genome];
            let mutations = &mut self.mutations[index_range.start..index_range.stop];
            ParentalGenomeMut {
                mutations,
                current_mutation_index: 0,
                genome,
            }
        } else {
            ParentalGenomeMut {
                mutations: &mut [],
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

fn compact_mutations(pop: &mut DiploidPopulation) {
    let new_indexes = {
        let mut temp = vec![];
        let mut i = 0_u32;
        temp.resize_with(pop.mutations.len(), || {
            let j = i;
            i += 1;
            j
        });
        let twon = 2 * (pop.individuals.len() as u32);
        let (mut left, _): (Vec<u32>, Vec<u32>) = temp.into_iter().partition(|&index| {
            let mcount = *unsafe { pop.mutation_counts.get_unchecked(index as usize) };
            mcount > 0 && mcount < twon
        });
        left.sort_unstable_by_key(|&m| {
            unsafe { pop.mutations.get_unchecked(m as usize) }.position()
        });
        left
    };
    let reindex = {
        let mut temp = vec![u32::MAX; pop.mutations.len()];
        for i in 0..new_indexes.len() {
            temp[new_indexes[i] as usize] = i as u32;
        }
        temp
    };

    // remap the genomes
    for m in &mut pop.genomes.mutations {
        *m = reindex[*m as usize];
    }

    let mut mutations = vec![];
    std::mem::swap(&mut mutations, &mut pop.mutations);
    // convert all mutations into options
    let mut mutations = mutations.into_iter().map(Some).collect::<Vec<_>>();
    let mut mutation_counts = vec![];
    for i in new_indexes {
        mutation_counts.push(pop.mutation_counts[i as usize]);
        pop.mutations.push(mutations[i as usize].take().unwrap());
    }
    std::mem::swap(&mut mutation_counts, &mut pop.mutation_counts);
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
    // NOTE: the Fn generic is for API flexibility
    // It is not always the case that we want to remove fixations:
    // * fixations affect phenotype (unlike popgen models of multiplicative fitness)
    // * multiple hits can make mutations "un-fixed" at a site
    fn mutation_recycling<F>(&self, f: F) -> Vec<usize>
    where
        F: Fn((usize, &u32)) -> Option<usize>,
    {
        self.mutation_counts
            .iter()
            .enumerate()
            .filter_map(f)
            .collect::<Vec<usize>>()
    }

    #[inline(never)]
    pub fn count_mutations(&mut self) {
        self.mutation_counts.fill(0);
        self.mutation_counts.resize(self.mutations.len(), 0);
        self.genomes
            .mutations
            .iter()
            // SAFETY: same issues as below
            .for_each(|m| *unsafe { self.mutation_counts.get_unchecked_mut(*m as usize) } += 1);
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

fn erase(mutation_counts: &[u32], twon: u32, mutations: &mut [u32]) -> usize {
    let mut rv = 0;

    // We only call this if we KNOW that the first mutation is a fixation,
    // so we can start to index from 1.
    for i in 1..mutations.len() {
        if *unsafe { mutation_counts.get_unchecked(mutations[i] as usize) } != twon {
            *unsafe { mutations.get_unchecked_mut(rv) } = mutations[i];
            rv += 1;
        }
    }
    rv
}

// This normally wouldn't be pub,
// but should be unit tested anyways.
#[inline(never)]
fn fixation_removal_check(mutation_counts: &[u32], twon: u32, input: &mut HaploidGenomes) -> bool {
    if mutation_counts.iter().any(|m| *m == twon) {
        let first_fixation = {
            let genome = input.get_genome(0);
            let index = genome
                .mutations
                .iter()
                .position(|&k| mutation_counts[k as usize] == twon)
                .unwrap();
            genome.mutations[index]
        };
        //let mut start = 0_usize;
        for i in 0..input.genome_spans.len() {
            let genome = input.get_genome_mut(i);
            let first_fixation_position = genome
                .mutations
                .iter()
                .position(|&m| m == first_fixation)
                .unwrap();
            let tpoint = erase(
                mutation_counts,
                twon,
                &mut genome.mutations[first_fixation_position..],
            );
            input.genome_spans[i].stop -=
                (genome.mutations.len() - tpoint - first_fixation_position);
            //output
            //    .mutations
            //    .extend_from_slice(&genome.mutations[..first_fixation_position]);
            //output.mutations.extend(
            //    genome
            //        .mutations
            //        .iter()
            //        .skip(first_fixation_position + 1)
            //        .filter(|&&m| {
            //            let count = *unsafe { mutation_counts.get_unchecked(m as usize) };
            //            count < twon
            //        }),
            //);
            // NOTE: this push can contradict what "add_range" is doing
            // in the impl of HaploidGenomes.
            //output.genome_spans.push(MutationRange {
            //    start,
            //    stop: output.mutations.len(),
            //});
            //start = output.mutations.len();
        }
        true
    } else {
        false
    }

    //if mutation_counts.iter().any(|m| *m == twon) {
    //    // NOTE: using this as an assert
    //    // eliminates all performance
    //    // benefits of the unsafe code below
    //    debug_assert!(output
    //        .mutations
    //        .iter()
    //        .all(|m| (*m as usize) < mutation_counts.len()));
    //    let x = output.mutations.len();
    //    // SAFETY: None in general. :(
    //    // In the specific case of this function not being public,
    //    // genomes only contain mutation keys generated
    //    // as they are added to the pop's mutation array.
    //    // As the same time, we make sure that mutation_counts'
    //    // length is also correct.
    //    //
    //    // This is a sticky issue:
    //    // 1. The perf gain of unsafe here is quite big.
    //    // 2. But it is gonna be hard to encapsulate things
    //    //    to the point where we feel generally good about saftety
    //    // 3. Therefore, it seems that the function itself should
    //    //    be marked unsafe in the absence of such encapsulation.
    //    output
    //        .mutations
    //        .retain(|m| *unsafe { mutation_counts.get_unchecked(*m as usize) } < twon);
    //    let delta = x - output.mutations.len();
    //    assert_eq!(delta % output.genome_spans.len(), 0);

    //    let delta_per_genome = delta / output.genome_spans.len();

    //    // NOTE: could be SIMD later
    //    output
    //        .genome_spans
    //        .iter_mut()
    //        .enumerate()
    //        .for_each(|(i, h)| {
    //            let c = h.start;
    //            h.start -= i * delta_per_genome;
    //            if i > 0 {
    //                assert!(c > h.start);
    //            }
    //            assert!(h.start < output.mutations.len());
    //            h.stop -= (i + 1) * delta_per_genome;
    //            assert!(h.stop >= h.start);
    //            assert!(
    //                h.stop <= output.mutations.len(),
    //                "{h:?} {}",
    //                output.mutations.len()
    //            );
    //        });
    //    true
    //} else {
    //    false
    //}
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
    todo!("implementing this nicely is all but impossible? the borrow checker hates what we are doing");
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
    let mut queue: Vec<usize> = vec![];
    for generation in 0..params.num_generations {
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
        
        // HACK to let us check performance
        pop.genomes.clear();
        pop.genomes
            .mutations
            .extend_from_slice(&offspring_genomes.mutations);
        pop.genomes
            .genome_spans
            .extend_from_slice(&offspring_genomes.genome_spans);
        offspring_genomes.clear();

        std::mem::swap(&mut pop.individuals, &mut offspring);
        offspring.clear();
        if generation % params.gcinterval == 0 {
            pop.count_mutations();

            if fixation_removal_check(
                &pop.mutation_counts,
                2 * params.num_individuals,
                &mut pop.genomes,
            ) {
                //std::mem::swap(&mut pop.genomes, &mut offspring_genomes);
                //offspring_genomes.clear();
            }
            queue = pop.mutation_recycling(|(index, &value)| {
                if value == 0 || value == 2 * params.num_individuals {
                    Some(index)
                } else {
                    None
                }
            });
        }
        if generation % 100 == 0 {
            compact_mutations(&mut pop);
            queue = vec![];
        }
    }
    Some(pop)
}
