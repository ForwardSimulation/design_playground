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

#[derive(Debug)]
struct Haplotypes {
    haplotypes: Vec<MutationRange>,
    mutations: Vec<usize>,
}

// FIXME: We are not using the type system
// well here.
// We need a way for a DiploidGenome
// to indicate that one of its elements
// is mutation-free.  Option<usize>
// seems the best here, but costs 2x the storage.
impl Default for Haplotypes {
    fn default() -> Self {
        let haplotypes = vec![MutationRange { start: 0, stop: 0 }];
        Self {
            haplotypes,
            mutations: vec![],
        }
    }
}

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
            let mut haplotypes = Haplotypes::default();

            // Now, everyone starts with a single "empty"
            // genome
            let individuals = vec![DiploidGenome::new(0, 0); size as usize];

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
//
// FIXME: code duplication...
fn make_offspring_genome(
    parent: DiploidGenome,
    parent_haplotypes: &Haplotypes,
    mutations: &Vec<Mutation>,
    new_mutations: Vec<usize>,
    breakpoints: &[Breakpoint],
    offspring_haplotypes: &mut Haplotypes,
) -> usize {
    assert!(breakpoints.is_empty(), "not dealing with recombination yet");
    let parent_range = parent_haplotypes.haplotypes[parent.first];
    if new_mutations.is_empty() {
        // TODO: "mendel" needs to happen here: 1/2 the time
        // we should use parent.second
        if parent_range.stop > parent_range.start {
            let parent_slice = &parent_haplotypes.mutations[parent_range.start..parent_range.stop];
            let start = offspring_haplotypes.mutations.len();
            offspring_haplotypes
                .mutations
                .extend_from_slice(parent_slice);
            let stop = offspring_haplotypes.mutations.len();
            let rv = offspring_haplotypes.haplotypes.len();
            offspring_haplotypes
                .haplotypes
                .push(MutationRange { start, stop });
            return rv;
        } else {
            return 0;
        }
    }

    if parent_range.stop > parent_range.start {
        let parent_slice = &parent_haplotypes.mutations[parent_range.start..parent_range.stop];
        let mut i = 0;
        let start = offspring_haplotypes.mutations.len();
        let nm = new_mutations.len();
        for m in new_mutations.iter() {
            let n = parent_slice[i..]
                .iter()
                .take_while(|mutation| mutations[**mutation].position() < mutations[*m].position())
                .inspect(|x| offspring_haplotypes.mutations.push(**x))
                .count();
            offspring_haplotypes.mutations.push(*m);
            i += n;
        }
        parent_slice[i..]
            .iter()
            .for_each(|m| offspring_haplotypes.mutations.push(*m));
        let stop = offspring_haplotypes.mutations.len();
        assert_eq!(
            stop - start,
            nm + parent_slice.len(),
            "{parent_slice:?} + {new_mutations:?} = {:?}",
            &offspring_haplotypes.mutations[start..stop]
        );
        let rv = offspring_haplotypes.haplotypes.len();
        offspring_haplotypes
            .haplotypes
            .push(MutationRange { start, stop });
        rv
    } else {
        if new_mutations.is_empty() {
            return 0;
        }
        let start = offspring_haplotypes.mutations.len();
        for m in new_mutations.into_iter() {
            offspring_haplotypes.mutations.push(m);
        }
        let stop = offspring_haplotypes.mutations.len();
        let rv = offspring_haplotypes.haplotypes.len();
        offspring_haplotypes
            .haplotypes
            .push(MutationRange { start, stop });
        rv
    }
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

    for generation in 0..params.num_generations {
        let mut offspring_haplotypes = Haplotypes::default();
        let mut offspring = vec![];
        let mut queue = pop.mutation_recycling();
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
            let first = make_offspring_genome(
                pop.individuals[parent1],
                &pop.haplotypes,
                &pop.mutations,
                mutations,
                &[], // no recombination for now...
                &mut offspring_haplotypes,
            );

            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut queue,
                &mut pop.mutations,
                &mut rng,
            );

            let second = make_offspring_genome(
                pop.individuals[parent2],
                &pop.haplotypes,
                &pop.mutations,
                mutations,
                &[], // no recombination for now...
                &mut offspring_haplotypes,
            );
            offspring.push(DiploidGenome { first, second });
        }
        pop.haplotypes = offspring_haplotypes;
        pop.individuals = offspring;
        println!("done with {generation}");
    }

    Some(pop)
}

#[test]
fn run_sim() {
    let params = SimParams {
        seed: 666,
        size: 1000,
        num_generations: 1000,
        mutation_rate: 1e-1,
    };
    // Empty genetic map == no recombination
    let builder = forrustts::genetics::GeneticMapBuilder::default();
    let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
    let _ = evolve_pop_with_haplotypes(params, genetic_map).unwrap();
}