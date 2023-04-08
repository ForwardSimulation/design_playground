use rand::prelude::Rng;
use rand::SeedableRng;

use forrustts::genetics::{Breakpoint, GeneticMap};
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
}

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

// NOTE: much of the behavior
// here should be associated fns
// of various types and/or other
// standalone fns.
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

            // ignore recombination and Mendel for now
            // and only pass on the 1st genome from
            // a parent + mutations
            let first = generate_offspring_genome(
                pop.individuals[parent1],
                &pop.haplotypes,
                &pop.mutations,
                mutations,
                &[], // no recombination for now...
                &mut offspring_haplotypes,
                &mut rng,
            );

            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut queue,
                &mut pop.mutations,
                &mut rng,
            );

            let second = generate_offspring_genome(
                pop.individuals[parent2],
                &pop.haplotypes,
                &pop.mutations,
                mutations,
                &[], // no recombination for now...
                &mut offspring_haplotypes,
                &mut rng,
            );
            offspring.push(DiploidGenome { first, second });
        }
        pop.haplotypes = offspring_haplotypes;
        pop.individuals = offspring;
        pop.count_mutations();
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
}
