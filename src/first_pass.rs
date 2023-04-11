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
use crate::common::SimParams;

#[derive(Debug, Default)]
struct Haplotypes {
    haplotypes: Vec<MutationRange>,
    mutations: Vec<usize>,
}

impl Haplotypes {
    fn get_genome(&self, genome: usize) -> ParentalGenome {
        if genome != usize::MAX {
            let index_range = self.haplotypes[genome];
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
            let rv = self.haplotypes.len();
            self.haplotypes.push(range);
            rv
        }
    }
}

// Implementation is specific to "diploid",
// so not an associated fn of Haplotypes
fn get_parental_genomes(
    haplotypes: &Haplotypes,
    parent: DiploidGenome,
) -> (ParentalGenome, ParentalGenome) {
    (
        haplotypes.get_genome(parent.first),
        haplotypes.get_genome(parent.second),
    )
}

/// When will the borrow checker hate this?
pub struct DiploidPopWithHaplotypes {
    haplotypes: Haplotypes,
    individuals: Vec<DiploidGenome>,
    mutations: Vec<Mutation>,
    mutation_counts: Vec<u32>,
}

impl DiploidPopWithHaplotypes {
    pub fn new(size: u32) -> Option<Self> {
        if size > 0 {
            let haplotypes = Haplotypes::default();

            // Now, everyone starts with a single "empty"
            // genome
            let individuals = vec![DiploidGenome::new(usize::MAX, usize::MAX); size as usize];

            Some(Self {
                haplotypes,
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
        self.haplotypes
            .mutations
            .iter()
            .for_each(|m| self.mutation_counts[*m] += 1);
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

// A proper implementation
// would be generic over "generating mutations"
#[inline(never)]
pub fn evolve_pop_with_haplotypes(
    params: SimParams,
    genetic_map: GeneticMap,
) -> Option<DiploidPopWithHaplotypes> {
    let params = params.validate()?;
    let mut pop = DiploidPopWithHaplotypes::new(params.num_individuals)?;

    let mut rng = rand::rngs::StdRng::seed_from_u64(params.seed);
    let parent_picker =
        rand::distributions::Uniform::<usize>::new(0, params.num_individuals as usize);
    let num_mutations = rand_distr::Poisson::<f64>::new(params.mutation_rate).ok()?;
    let position_generator = rand::distributions::Uniform::<Position>::new(
        Position::new_valid(0),
        Position::new_valid(1000000),
    );
    let u01 = rand::distributions::Uniform::new(0., 1.);

    let mut genetic_map = genetic_map;

    //let mut parent_haplotype_map = vec![];
    let mut offspring_haplotypes = Haplotypes::default();
    let mut offspring = vec![];
    for generation in 0..params.num_generations {
        offspring_haplotypes.mutations.reserve(1000);
        let mut queue = pop.mutation_recycling();
        for _ in 0..params.num_individuals {
            // Pick two parents
            let parent1 = rng.sample(parent_picker);
            let parent2 = rng.sample(parent_picker);

            // Mutations for offspring genome 1
            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut queue,
                &mut pop.mutations,
                &mut rng,
            );

            genetic_map.generate_breakpoints(&mut rng);
            pop.mutation_counts
                .resize(pop.mutation_counts.len() + mutations.len(), 0);

            let genomes = get_parental_genomes(&pop.haplotypes, pop.individuals[parent1]);
            let genomes = if rng.sample(u01) < 0.5 {
                genomes
            } else {
                (genomes.1, genomes.0)
            };
            let range = generate_offspring_genome(
                genomes,
                &pop.mutations,
                &pop.mutation_counts,
                2 * params.num_individuals,
                mutations,
                genetic_map.breakpoints(),
                &mut offspring_haplotypes.mutations,
            );
            let first = offspring_haplotypes.add_range(range);

            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut queue,
                &mut pop.mutations,
                &mut rng,
            );

            genetic_map.generate_breakpoints(&mut rng);
            pop.mutation_counts
                .resize(pop.mutation_counts.len() + mutations.len(), 0);

            let genomes = get_parental_genomes(&pop.haplotypes, pop.individuals[parent2]);
            let genomes = if rng.sample(u01) < 0.5 {
                genomes
            } else {
                (genomes.1, genomes.0)
            };
            let range = generate_offspring_genome(
                genomes,
                &pop.mutations,
                &pop.mutation_counts,
                2 * params.num_individuals,
                mutations,
                genetic_map.breakpoints(),
                &mut offspring_haplotypes.mutations,
            );
            let second = offspring_haplotypes.add_range(range);
            offspring.push(DiploidGenome { first, second });
        }
        std::mem::swap(&mut pop.haplotypes, &mut offspring_haplotypes);
        offspring_haplotypes.haplotypes.clear();
        offspring_haplotypes.mutations.clear();

        std::mem::swap(&mut pop.individuals, &mut offspring);
        offspring.clear();
        pop.count_mutations();
        // for h in &pop.haplotypes.haplotypes {
        //     assert!(pop.haplotypes.mutations[h.start..h.stop]
        //         .windows(2)
        //         .all(|w| pop.mutations[w[0]].position() <= pop.mutations[w[1]].position()));
        // }
        // println!(
        //     "{} {}/{}",
        //     params.mutation_rate,
        //     parent_haplotype_map
        //         .iter()
        //         .filter(|i| **i != usize::MAX)
        //         .count(),
        //     pop.haplotypes.haplotypes.len(),
        // );
        //println!("done with {generation}, {}", pop.mutations.len());
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
            let _ = evolve_pop_with_haplotypes(params, genetic_map).unwrap();
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
            // Empty genetic map == no recombination
            let builder = forrustts::genetics::GeneticMapBuilder::default().extend_poisson(&poisson);

            let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
            let _ = evolve_pop_with_haplotypes(params, genetic_map).unwrap();
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
