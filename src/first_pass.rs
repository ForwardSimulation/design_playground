use clap::Parser;
use rand::prelude::Rng;
use rand::SeedableRng;

use forrustts::genetics::{Breakpoint, GenerateBreakpoints, GeneticMap};
use forrustts::prelude::*;

// We need a type with a more complex
// layout than a simple position.
// It is > 32 bits just to try
// to annoy memory accesses.
#[derive(Clone, Debug)]
pub struct Mutation {
    position: Position,
    effect_sizes: Vec<f64>,
    origin_time: Time,
}

impl Mutation {
    pub fn new(position: Position, effect_sizes: Vec<f64>, origin_time: Time) -> Self {
        Self {
            position,
            effect_sizes,
            origin_time,
        }
    }

    pub fn position(&self) -> Position {
        self.position
    }

    pub fn effect_sizes(&self) -> &[f64] {
        &self.effect_sizes
    }
}

#[derive(Copy, Clone, Debug)]
struct MutationRange {
    start: usize,
    stop: usize,
}

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

#[derive(Copy, Clone, Debug)]
struct ParentalGenome<'a> {
    mutations: &'a [usize],
    current_mutation_index: usize,
    genome: usize,
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

// NOTE: we use usize::MAX to indicate a
// "no genome" state. Production
// code would do better.
#[derive(Copy, Clone, Debug)]
struct DiploidGenome {
    first: usize,
    second: usize,
}

impl DiploidGenome {
    pub fn new(first: usize, second: usize) -> Self {
        Self { first, second }
    }
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

#[derive(Parser, Debug)]
pub struct SimParams {
    #[arg(short, long)]
    pub seed: u64,
    #[arg(short, long)]
    pub size: u32,
    #[arg(short, long)]
    pub num_generations: u32,
    #[arg(short, long)]
    pub mutation_rate: f64,
    #[arg(short, long)]
    pub recrate: f64,
}

// NOTE: we can be smarter than this
// by using a struct that can populate
// a Vec and output a slice, thus avoiding
// one allocation per offspring.
#[inline(never)]
fn generate_mutations(
    generation: u32,
    num_mutations: rand_distr::Poisson<f64>,
    position_generator: rand::distributions::Uniform<Position>,
    queue: &mut Vec<usize>,
    mutations: &mut Vec<Mutation>,
    rng: &mut rand::rngs::StdRng,
) -> Vec<usize> {
    let mut rv = vec![];

    // Dangerous if mean is very large,
    // but of little practical consequence
    let nmuts = rng.sample(num_mutations) as u32;

    for _ in 0..nmuts {
        let pos = rng.sample(position_generator);
        // Insufficiently generic -- we are really
        // just making Positions here, and there is
        // no information about how to generate
        // other fields.
        let mutation = Mutation::new(pos, vec![0.0], Time::from(generation as i64));
        match queue.pop() {
            Some(index) => {
                mutations[index] = mutation;
                rv.push(index);
            }
            None => {
                rv.push(mutations.len());
                mutations.push(mutation);
            }
        };
    }

    // sort our new mutations by (increasing) position
    rv.sort_by(|i, j| mutations[*i].position().cmp(&mutations[*j].position()));

    rv
}

#[inline(never)]
fn merge_mutations(
    mutations: &[Mutation],
    new_mutations: &[usize],
    offspring_haplotypes: &mut Vec<usize>,
    current_genome: &mut ParentalGenome,
) {
    for m in new_mutations.iter() {
        let n = current_genome.mutations[current_genome.current_mutation_index..]
            .iter()
            .take_while(|mutation| mutations[**mutation].position() < mutations[*m].position())
            .inspect(|x| offspring_haplotypes.push(**x))
            .count();
        offspring_haplotypes.push(*m);
        current_genome.current_mutation_index += n;
    }
    offspring_haplotypes
        .extend_from_slice(&current_genome.mutations[current_genome.current_mutation_index..]);
}

// This has had a lot of refactoring and still
// is hard to test in isolation (see tests below).
#[inline(never)]
fn generate_offspring_genome(
    genomes: (ParentalGenome, ParentalGenome),
    mutations: &[Mutation],
    new_mutations: Vec<usize>,
    breakpoints: &[Breakpoint],
    offspring_mutations: &mut Vec<usize>,
) -> MutationRange {
    let (mut current_genome, mut other_genome) = genomes;
    let start = offspring_mutations.len();
    if breakpoints.is_empty() {
        merge_mutations(
            mutations,
            &new_mutations,
            offspring_mutations,
            &mut current_genome,
        );
    } else {
        let mut mut_index = 0_usize;
        for b in breakpoints {
            let bpos = match b {
                // NOTE: forrustts needs to handle this
                // comparison with trat impls
                forrustts::genetics::Breakpoint::Crossover(pos) => *pos,
                forrustts::genetics::Breakpoint::IndependentAssortment(pos) => *pos,
                _ => unimplemented!("unhandled Breakpoint variant"),
            };
            mut_index += new_mutations[mut_index..]
                .iter()
                .take_while(|k| mutations[**k].position() < bpos)
                .inspect(|k| {
                    current_genome.current_mutation_index += current_genome.mutations
                        [current_genome.current_mutation_index..]
                        .iter()
                        .take_while(|gk| mutations[**gk].position() < mutations[**k].position())
                        .inspect(|gk| offspring_mutations.push(**gk))
                        .count();
                    offspring_mutations.push(**k);
                })
                .count();
            current_genome.current_mutation_index += current_genome.mutations
                [current_genome.current_mutation_index..]
                .iter()
                .take_while(|gk| mutations[**gk].position() < bpos)
                .inspect(|gk| offspring_mutations.push(**gk))
                .count();

            // Advance other genome
            other_genome.current_mutation_index += other_genome.mutations
                [other_genome.current_mutation_index..]
                .iter()
                .take_while(|gk| mutations[**gk].position() < bpos)
                .count();

            std::mem::swap(&mut current_genome, &mut other_genome);
        }
        merge_mutations(
            mutations,
            &new_mutations[mut_index..],
            offspring_mutations,
            &mut current_genome,
        )
    }
    let stop = offspring_mutations.len();
    MutationRange { start, stop }
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
    let mut pop = DiploidPopWithHaplotypes::new(params.size)?;

    let mut rng = rand::rngs::StdRng::seed_from_u64(params.seed);
    let parent_picker = rand::distributions::Uniform::<usize>::new(0, params.size as usize);
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
        for _ in 0..params.size {
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

            let genomes = get_parental_genomes(&pop.haplotypes, pop.individuals[parent1]);
            let genomes = if rng.sample(u01) < 0.5 {
                genomes
            } else {
                (genomes.1, genomes.0)
            };
            let range = generate_offspring_genome(
                genomes,
                &pop.mutations,
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

            let genomes = get_parental_genomes(&pop.haplotypes, pop.individuals[parent2]);
            let genomes = if rng.sample(u01) < 0.5 {
                genomes
            } else {
                (genomes.1, genomes.0)
            };
            let range = generate_offspring_genome(
                genomes,
                &pop.mutations,
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
                size: 100,
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
                size: 100,
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
mod test_create_offspring_genome {
    use super::*;
    use proptest::prelude::*;

    fn setup(
        seed: u64,
        nmuts1: usize,
        nmuts2: usize,
        num_new_mutations: usize,
        nbreakpoints: usize,
    ) -> (
        Vec<Mutation>,
        Vec<usize>,
        Vec<usize>,
        Vec<forrustts::genetics::Breakpoint>,
    ) {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let start = forrustts::Position::new_valid(0);
        let stop = forrustts::Position::new_valid(100000);
        let makepos = rand::distributions::Uniform::new(start, stop);

        let mut mutations = vec![];
        let mut haploid_genomes = vec![];
        for _ in 0..(nmuts1 + nmuts2) {
            let position = rng.sample(makepos);
            let m = Mutation {
                position,
                effect_sizes: vec![0.1],
                origin_time: 0.into(),
            };
            haploid_genomes.push(mutations.len());
            mutations.push(m);
        }

        let mut new_mutations = vec![];
        for _ in 0..(num_new_mutations) {
            let position = rng.sample(makepos);
            let m = Mutation {
                position,
                effect_sizes: vec![0.1],
                origin_time: 0.into(),
            };
            new_mutations.push(mutations.len());
            mutations.push(m);
        }

        let mut breakpoints = vec![];
        for _ in 0..nbreakpoints {
            let position = rng.sample(makepos);
            breakpoints.push(Breakpoint::Crossover(position));
        }

        new_mutations.sort_by(|i, j| mutations[*i].position().cmp(&mutations[*j].position()));
        breakpoints.sort();
        haploid_genomes[0..nmuts1]
            .sort_by(|i, j| mutations[*i].position().cmp(&mutations[*j].position()));
        haploid_genomes[nmuts1..]
            .sort_by(|i, j| mutations[*i].position().cmp(&mutations[*j].position()));
        (mutations, haploid_genomes, new_mutations, breakpoints)
    }

    fn setup_parents(
        nmuts1: usize,
        nmuts2: usize,
        haploid_genomes: &[usize],
    ) -> (ParentalGenome, ParentalGenome) {
        let parent1_genome = if nmuts1 > 0 {
            ParentalGenome {
                mutations: &haploid_genomes[0..nmuts1],
                current_mutation_index: 0,
                genome: usize::MAX,
            }
        } else {
            ParentalGenome {
                mutations: &[],
                current_mutation_index: 0,
                genome: usize::MAX,
            }
        };
        let parent2_genome = if nmuts1 > 0 {
            ParentalGenome {
                mutations: &haploid_genomes[nmuts1..(nmuts1 + nmuts2)],
                current_mutation_index: 0,
                genome: usize::MAX,
            }
        } else {
            ParentalGenome {
                mutations: &[],
                current_mutation_index: 0,
                genome: usize::MAX,
            }
        };
        (parent1_genome, parent2_genome)
    }

    fn naive(
        parents: (ParentalGenome, ParentalGenome),
        mutations: &[Mutation],
        breakpoints: &[Breakpoint],
        new_mutations: &[usize],
    ) -> Vec<usize> {
        let mut kept_breakpoints = vec![];
        for b in breakpoints.iter() {
            let x = breakpoints.iter().filter(|&i| b == i).count();
            if x % 2 != 0 {
                kept_breakpoints.push(*b);
            }
        }
        let (mut parent1_genome, mut parent2_genome) = parents;
        let mut lastpos = Position::new_valid(0);
        let mut output = vec![];

        for b in kept_breakpoints.iter() {
            let pos = match b {
                Breakpoint::Crossover(x) => x,
                Breakpoint::IndependentAssortment(x) => x,
                _ => unimplemented!("bad"),
            };
            parent1_genome
                .mutations
                .iter()
                .filter(|&k| mutations[*k].position() >= lastpos && mutations[*k].position() < *pos)
                .for_each(|a| output.push(*a));
            lastpos = *pos;
            std::mem::swap(&mut parent1_genome, &mut parent2_genome);
        }
        parent1_genome
            .mutations
            .iter()
            .filter(|&k| mutations[*k].position() >= lastpos)
            .for_each(|a| output.push(*a));
        for m in new_mutations {
            output.push(*m);
        }
        output.sort_by(|i, j| mutations[*i].position().cmp(&mutations[*j].position()));
        output
    }

    fn keys_to_positions(mutations: &[Mutation], keys: &[usize]) -> Vec<i64> {
        let mut rv = vec![];

        for k in keys {
            rv.push(mutations[*k].position().into())
        }

        rv
    }

    fn run(seed: u64, nmuts1: usize, nmuts2: usize, num_new_mutations: usize, nbreakpoints: usize) {
        let (mutations, haploid_genomes, new_mutations, breakpoints) =
            setup(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
        let (parent1_genome, parent2_genome) = setup_parents(nmuts1, nmuts2, &haploid_genomes);
        assert!(parent1_genome
            .mutations
            .windows(2)
            .all(|w| mutations[w[0]].position() <= mutations[w[1]].position()),);
        assert!(parent2_genome
            .mutations
            .windows(2)
            .all(|w| mutations[w[0]].position() <= mutations[w[1]].position()),);
        let mut offspring_genomes = Vec::<usize>::new();
        let range = generate_offspring_genome(
            (parent1_genome, parent2_genome),
            &mutations,
            new_mutations.clone(),
            &breakpoints,
            &mut offspring_genomes,
        );
        let naive_output = naive(
            (parent1_genome, parent2_genome),
            &mutations,
            &breakpoints,
            &new_mutations,
        );
        assert!(naive_output
            .windows(2)
            .all(|w| mutations[w[0]].position() <= mutations[w[1]].position()),);
        assert_eq!(range.stop, offspring_genomes.len());
        assert_eq!(
            naive_output.len(),
            offspring_genomes.len(),
            "[{:?}, {:?}] + {:?} & {:?} = {:?} and {:?}",
            keys_to_positions(&mutations, parent1_genome.mutations),
            keys_to_positions(&mutations, parent2_genome.mutations),
            keys_to_positions(&mutations, &new_mutations),
            breakpoints,
            keys_to_positions(&mutations, &naive_output),
            keys_to_positions(&mutations, &offspring_genomes),
        );
        for i in &offspring_genomes {
            assert_eq!(naive_output.iter().filter(|&j| j == i).count(), 1);
        }
    }

    proptest! {
        #[test]
        fn test_with_no_breakpoints(seed in 0..u64::MAX,
                                    nmuts1 in 0..=50_usize,
                                    nmuts2 in 0..=50_usize,
                                    num_new_mutations in 0..50_usize,
                                    )
        {
            run(seed, nmuts1, nmuts2, num_new_mutations, 0);
        }

        #[test]
        fn test_with_breakpoints(seed in 0..u64::MAX,
                                    nmuts1 in 0..=50_usize,
                                    nmuts2 in 0..=50_usize,
                                    num_new_mutations in 0..50_usize,
                                    nbreakpoints in 0..25_usize
                                    )
        {
            run(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
        }
    }

    #[test]
    fn proptest_regression_0() {
        run(0, 0, 0, 1, 0);
    }

    #[test]
    fn proptest_regression_1() {
        let seed = 3958894411756207915_u64;
        let nmuts1 = 28;
        let nmuts2 = 32;
        let num_new_mutations = 30;
        run(seed, nmuts1, nmuts2, num_new_mutations, 0);
    }

    #[test]
    fn proptest_regression_2() {
        let seed = 0;
        let nmuts1 = 1;
        let nmuts2 = 0;
        let num_new_mutations = 8;
        let nbreakpoints = 1;
        run(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
    }

    #[test]
    fn proptest_regression_3() {
        let seed = 7732215647961344_u64;
        let nmuts1 = 0;
        let nmuts2 = 17;
        let num_new_mutations = 14;
        let nbreakpoints = 1;
        run(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
    }

    #[test]
    fn proptest_regression_4() {
        let seed = 2727126891865037407_u64;
        let nmuts1 = 20;
        let nmuts2 = 11;
        let num_new_mutations = 4;
        let nbreakpoints = 1;
        run(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
    }

    #[test]
    fn proptest_regression_5() {
        let seed = 7956035;
        let nmuts1 = 3;
        let nmuts2 = 11;
        let num_new_mutations = 0;
        let nbreakpoints = 1;
        run(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
    }

    #[test]
    fn proptest_regression_6() {
        let seed = 0;
        let nmuts1 = 1;
        let nmuts2 = 10;
        let num_new_mutations = 0;
        let nbreakpoints = 4;
        run(seed, nmuts1, nmuts2, num_new_mutations, nbreakpoints);
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
