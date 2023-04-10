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
    pub fn mutation_recycling(&self) -> Vec<usize> {
        self.mutation_counts
            .iter()
            .enumerate()
            .filter_map(|(index, value)| if value == &0 { Some(index) } else { None })
            .collect::<Vec<usize>>()
    }

    pub fn count_mutations(&mut self) {
        self.mutation_counts.fill(0);
        self.mutation_counts.resize(self.mutations.len(), 0);
        self.haplotypes
            .mutations
            .iter()
            .for_each(|m| self.mutation_counts[*m] += 1);
    }
}

pub struct SimParams {
    pub seed: u64,
    pub size: u32,
    pub num_generations: u32,
    pub mutation_rate: f64,
}

// NOTE: we can be smarter than this
// by using a struct that can populate
// a Vec and output a slice, thus avoiding
// one allocation per offspring.
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

fn generate_offspring_genome(
    parent: DiploidGenome,
    parent_haplotypes: &Haplotypes,
    mutations: &[Mutation],
    new_mutations: Vec<usize>,
    breakpoints: &[Breakpoint],
    offspring_haplotypes: &mut Haplotypes,
    rng: &mut rand::rngs::StdRng,
) -> usize {
    let u01 = rand::distributions::Uniform::new(0., 1.);
    let (mut first_genome, mut second_genome) = get_parental_genomes(parent_haplotypes, parent);
    if rng.sample(u01) < 0.5 {
        std::mem::swap(&mut first_genome, &mut second_genome);
    }
    let mut rv = usize::MAX;
    let start = offspring_haplotypes.mutations.len();
    let nm = new_mutations.len();
    for m in new_mutations.iter() {
        let n = first_genome.mutations[first_genome.current_mutation_index..]
            .iter()
            .take_while(|mutation| mutations[**mutation].position() < mutations[*m].position())
            .inspect(|x| offspring_haplotypes.mutations.push(**x))
            .count();
        offspring_haplotypes.mutations.push(*m);
        first_genome.current_mutation_index += n;
    }
    first_genome.mutations[first_genome.current_mutation_index..]
        .iter()
        .for_each(|m| offspring_haplotypes.mutations.push(*m));
    let stop = offspring_haplotypes.mutations.len();
    if stop > start {
        rv = offspring_haplotypes.haplotypes.len();
        offspring_haplotypes
            .haplotypes
            .push(MutationRange { start, stop });
    }
    assert_eq!(
        stop - start,
        nm + first_genome.mutations.len(),
        "{:?} + {new_mutations:?} = {:?}",
        first_genome.mutations,
        &offspring_haplotypes.mutations[start..stop]
    );
    rv
}

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
    current_genome.mutations[current_genome.current_mutation_index..]
        .iter()
        .for_each(|m| offspring_haplotypes.push(*m));
}

// Too much coupling makes all of this untestable.
fn generate_offspring_genome_test(
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
            println!("mut index = {}", mut_index);
            mut_index += new_mutations[mut_index..]
                .iter()
                .take_while(|k| match b {
                    // NOTE: forrustts needs to handle this
                    // comparison with trat impls
                    forrustts::genetics::Breakpoint::Crossover(pos) => {
                        mutations[**k].position() < *pos
                    }
                    forrustts::genetics::Breakpoint::IndependentAssortment(pos) => {
                        mutations[**k].position() < *pos
                    }
                    _ => unimplemented!("unhandled Breakpoint variant"),
                })
                .inspect(|k| {
                    current_genome.current_mutation_index += current_genome
                        .mutations
                        .iter()
                        .take_while(|gk| mutations[**gk].position() < mutations[**k].position())
                        .inspect(|gk| offspring_mutations.push(**gk))
                        .count();
                    offspring_mutations.push(**k);
                })
                .count();
            println!("mut index = {}", mut_index);
            println!(
                "{} {}, {} {}",
                current_genome.mutations.len(),
                current_genome.current_mutation_index,
                other_genome.mutations.len(),
                other_genome.current_mutation_index
            );
            //current_genome.current_mutation_index += current_genome.mutations
            //    [current_genome.current_mutation_index..]
            //    .iter()
            //    .take_while(|gk| match b {
            //        // NOTE: forrustts needs to handle this
            //        // comparison with trat impls
            //        forrustts::genetics::Breakpoint::Crossover(pos) => {
            //            mutations[**gk].position() < *pos
            //        }
            //        forrustts::genetics::Breakpoint::IndependentAssortment(pos) => {
            //            mutations[**gk].position() < *pos
            //        }
            //        _ => unimplemented!("unhandled Breakpoint variant"),
            //    })
            //    .inspect(|gk| offspring_mutations.push(**gk))
            //    .count();

            // Advance other genome
            other_genome.current_mutation_index += other_genome.mutations
                [other_genome.current_mutation_index..]
                .iter()
                .take_while(|gk| match b {
                    // NOTE: forrustts needs to handle this
                    // comparison with trat impls
                    forrustts::genetics::Breakpoint::Crossover(pos) => {
                        mutations[**gk].position() < *pos
                    }
                    forrustts::genetics::Breakpoint::IndependentAssortment(pos) => {
                        mutations[**gk].position() < *pos
                    }
                    _ => unimplemented!("unhandled Breakpoint variant"),
                })
                .count();

            std::mem::swap(&mut current_genome, &mut other_genome);
        }
        println!("mut index = {}", mut_index);
        println!(
            "{} {}, {} {}",
            current_genome.mutations.len(),
            current_genome.current_mutation_index,
            other_genome.mutations.len(),
            other_genome.current_mutation_index
        );
        //merge_mutations(
        //    mutations,
        //    &new_mutations[mut_index..],
        //    offspring_mutations,
        //    &mut current_genome,
        //)
    }
    let stop = offspring_mutations.len();
    MutationRange { start, stop }
    //if stop > start {
    //    rv = offspring_mutations.len();
    //    offspring_mutations
    //        .haplotypes
    //        .push(MutationRange { start, stop });
    //}
    //assert_eq!(
    //    stop - start,
    //    nm + current_genome.mutations.len(),
    //    "{:?} + {new_mutations:?} = {:?}",
    //    current_genome.mutations,
    //    &offspring_haplotypes.mutations[start..stop]
    //);
    //rv
}

fn generate_offspring_genome2(
    parent: DiploidGenome,
    parent_haplotypes: &Haplotypes,
    mutations: &[Mutation],
    new_mutations: Vec<usize>,
    breakpoints: &[Breakpoint],
    offspring_haplotypes: &mut Haplotypes,
    parent_haplotype_map: &mut [usize],
    rng: &mut rand::rngs::StdRng,
) -> usize {
    let u01 = rand::distributions::Uniform::new(0., 1.);
    let (mut first_genome, mut second_genome) = get_parental_genomes(parent_haplotypes, parent);
    if rng.sample(u01) < 0.5 {
        std::mem::swap(&mut first_genome, &mut second_genome);
    }
    let mut rv = usize::MAX;
    if breakpoints.is_empty() && new_mutations.is_empty() {
        // Then we have aready processed this genome
        if parent_haplotype_map[first_genome.genome] != usize::MAX {
            return parent_haplotype_map[first_genome.genome];
        }
    }
    let start = offspring_haplotypes.mutations.len();
    let nm = new_mutations.len();
    for m in new_mutations.iter() {
        let n = first_genome.mutations[first_genome.current_mutation_index..]
            .iter()
            .take_while(|mutation| mutations[**mutation].position() < mutations[*m].position())
            .inspect(|x| offspring_haplotypes.mutations.push(**x))
            .count();
        offspring_haplotypes.mutations.push(*m);
        first_genome.current_mutation_index += n;
    }
    first_genome.mutations[first_genome.current_mutation_index..]
        .iter()
        .for_each(|m| offspring_haplotypes.mutations.push(*m));
    let stop = offspring_haplotypes.mutations.len();
    if stop > start {
        rv = offspring_haplotypes.haplotypes.len();
        offspring_haplotypes
            .haplotypes
            .push(MutationRange { start, stop });
    }
    if nm == 0 && breakpoints.is_empty() {
        assert_eq!(parent_haplotype_map[first_genome.genome], usize::MAX);
        parent_haplotype_map[first_genome.genome] = rv;
    }
    assert_eq!(
        stop - start,
        nm + first_genome.mutations.len(),
        "{:?} + {new_mutations:?} = {:?}",
        first_genome.mutations,
        &offspring_haplotypes.mutations[start..stop]
    );
    rv
}

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
    for generation in 0..params.num_generations {
        let mut offspring_haplotypes = Haplotypes::default();
        let mut offspring = vec![];
        let mut queue = pop.mutation_recycling();
        //parent_haplotype_map.fill(usize::MAX);
        //parent_haplotype_map.resize(pop.haplotypes.haplotypes.len(), usize::MAX);
        for birth in 0..params.size {
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
            let range = generate_offspring_genome_test(
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
            let range = generate_offspring_genome_test(
                genomes,
                &pop.mutations,
                mutations,
                genetic_map.breakpoints(),
                &mut offspring_haplotypes.mutations,
            );
            let second = offspring_haplotypes.add_range(range);
            offspring.push(DiploidGenome { first, second });
        }
        pop.haplotypes = offspring_haplotypes;
        pop.individuals = offspring;
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
        fn run_sim_no_recombination(seed in 0..u64::MAX) {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            let make_mutrate = rand_distr::Exp::new(1.0).unwrap();
            let mutation_rate = rng.sample(make_mutrate);
            let params = SimParams {
                seed,
                size: 100,
                num_generations: 100,
                mutation_rate,
            };
            // Empty genetic map == no recombination
            let builder = forrustts::genetics::GeneticMapBuilder::default();
            let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
            let _ = evolve_pop_with_haplotypes(params, genetic_map).unwrap();
        }
    }

    proptest! {
        #[test]
        fn run_sim_with_recombination(seed in 0..u64::MAX) {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            let make_mutrate = rand_distr::Exp::new(1.0).unwrap();
            let mutation_rate = rng.sample(make_mutrate);
            let params = SimParams {
                seed,
                size: 100,
                num_generations: 100,
                mutation_rate,
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
    use proptest::{prelude::*, option::of};

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
            if x % 2 == 0 {
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
        let range = generate_offspring_genome_test(
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
            parent1_genome.mutations,
            parent2_genome.mutations,
            new_mutations,
            breakpoints,
            naive_output,
            offspring_genomes,
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
                                    nbreakpoints in 0..10_usize
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
