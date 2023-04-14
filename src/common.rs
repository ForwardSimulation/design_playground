use clap::Parser;
use rand::prelude::Rng;

use forrustts::genetics::Breakpoint;
use forrustts::prelude::*;

#[derive(Parser, Debug)]
pub struct SimParams {
    #[arg(short, long)]
    pub seed: u64,
    #[arg(long = "popsize", short = 'p')]
    pub num_individuals: u32,
    #[arg(short, long = "ngens")]
    pub num_generations: u32,
    #[arg(short, long = "mu")]
    pub mutation_rate: f64,
    #[arg(short, long = "r")]
    pub recrate: f64,
}

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
pub struct MutationRange {
    pub start: usize,
    pub stop: usize,
}

#[derive(Copy, Clone, Debug)]
pub struct ParentalGenome<'a> {
    pub mutations: &'a [u32],
    pub current_mutation_index: usize,
    pub genome: usize,
}

// NOTE: we use usize::MAX to indicate a
// "no genome" state. Production
// code would do better.
#[derive(Copy, Clone, Debug)]
pub struct DiploidGenome {
    pub first: usize,
    pub second: usize,
}

impl DiploidGenome {
    pub fn new(first: usize, second: usize) -> Self {
        Self { first, second }
    }
}

// NOTE: we can be smarter than this
// by using a struct that can populate
// a Vec and output a slice, thus avoiding
// one allocation per offspring.
#[inline(never)]
pub fn generate_mutations(
    generation: u32,
    num_mutations: rand_distr::Poisson<f64>,
    position_generator: rand::distributions::Uniform<Position>,
    queue: &mut Vec<usize>,
    mutations: &mut Vec<Mutation>,
    rng: &mut rand::rngs::StdRng,
) -> Vec<u32> {
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
                rv.push(index as u32);
            }
            None => {
                rv.push(mutations.len() as u32);
                mutations.push(mutation);
            }
        };
    }

    // sort our new mutations by (increasing) position
    rv.sort_by(|i, j| {
        mutations[*i as usize]
            .position()
            .cmp(&mutations[*j as usize].position())
    });

    rv
}

fn get_partition_index(genome: &[u32], mutations: &[Mutation], position: Position) -> usize {
    genome.partition_point(|k| mutations[*k as usize].position() < position)
}

fn extend_from_slice(
    genome: &[u32],
    mutations: &[Mutation],
    position: Position,
    offspring_haplotypes: &mut Vec<u32>,
) -> usize {
    let n = get_partition_index(genome, mutations, position);
    offspring_haplotypes.extend_from_slice(&genome[..(n as usize)]);
    n
}

fn update_genome(
    mutations: &[Mutation],
    position: Position,
    parent_genome: &mut ParentalGenome,
    offspring_haplotypes: &mut Vec<u32>,
) {
    parent_genome.current_mutation_index += extend_from_slice(
        &parent_genome.mutations[parent_genome.current_mutation_index..],
        mutations,
        position,
        offspring_haplotypes,
    )
}

#[inline(never)]
fn merge_mutations(
    mutations: &[Mutation],
    new_mutations: &[u32],
    offspring_haplotypes: &mut Vec<u32>,
    current_genome: &mut ParentalGenome,
) {
    for m in new_mutations.iter() {
        let mpos = mutations[*m as usize].position();
        update_genome(mutations, mpos, current_genome, offspring_haplotypes);
        offspring_haplotypes.push(*m);
    }
    offspring_haplotypes
        .extend_from_slice(&current_genome.mutations[current_genome.current_mutation_index..]);
}

// This has had a lot of refactoring and still
// is hard to test in isolation (see tests below).
#[inline(never)]
pub fn generate_offspring_genome(
    genomes: (ParentalGenome, ParentalGenome),
    mutations: &[Mutation],
    new_mutations: Vec<u32>,
    breakpoints: &[Breakpoint],
    offspring_mutations: &mut Vec<u32>,
) -> MutationRange {
    let (mut current_genome, mut other_genome) = genomes;
    let start = offspring_mutations.len();
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
            .take_while(|k| mutations[**k as usize].position() < bpos)
            .inspect(|k| {
                let mpos = mutations[**k as usize].position();
                update_genome(mutations, mpos, &mut current_genome, offspring_mutations);
                offspring_mutations.push(**k);
            })
            .count();
        update_genome(mutations, bpos, &mut current_genome, offspring_mutations);

        // Advance other genome
        other_genome.current_mutation_index += get_partition_index(
            &other_genome.mutations[other_genome.current_mutation_index..],
            mutations,
            bpos,
        );

        std::mem::swap(&mut current_genome, &mut other_genome);
    }
    merge_mutations(
        mutations,
        &new_mutations[mut_index..],
        offspring_mutations,
        &mut current_genome,
    );
    let stop = offspring_mutations.len();
    MutationRange { start, stop }
}

#[inline(never)]
pub fn set_fixation_counts_to_zero(twon: u32, mutation_counts: &mut [u32]) {
    mutation_counts
        .iter_mut()
        .filter(|m| **m == twon)
        .for_each(|m| *m = 0);
}

#[cfg(test)]
mod test_create_offspring_genome {
    use super::*;
    use proptest::prelude::*;
    use rand::SeedableRng;

    fn setup(
        seed: u64,
        nmuts1: usize,
        nmuts2: usize,
        num_new_mutations: usize,
        nbreakpoints: usize,
    ) -> (
        Vec<Mutation>,
        Vec<u32>,
        Vec<u32>,
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
            let m = Mutation::new(position, vec![0.1], 0.into());
            haploid_genomes.push(mutations.len() as u32);
            mutations.push(m);
        }

        let mut new_mutations = vec![];
        for _ in 0..(num_new_mutations) {
            let position = rng.sample(makepos);
            let m = Mutation::new(position, vec![0.1], 0.into());
            new_mutations.push(mutations.len() as u32);
            mutations.push(m);
        }

        let mut breakpoints = vec![];
        for _ in 0..nbreakpoints {
            let position = rng.sample(makepos);
            breakpoints.push(Breakpoint::Crossover(position));
        }

        new_mutations.sort_by(|i, j| {
            mutations[*i as usize]
                .position()
                .cmp(&mutations[*j as usize].position())
        });
        breakpoints.sort();
        haploid_genomes[0..nmuts1].sort_by(|i, j| {
            mutations[*i as usize]
                .position()
                .cmp(&mutations[*j as usize].position())
        });
        haploid_genomes[nmuts1..].sort_by(|i, j| {
            mutations[*i as usize]
                .position()
                .cmp(&mutations[*j as usize].position())
        });
        (mutations, haploid_genomes, new_mutations, breakpoints)
    }

    fn setup_parents(
        nmuts1: usize,
        nmuts2: usize,
        haploid_genomes: &[u32],
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
        new_mutations: &[u32],
    ) -> Vec<u32> {
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
                .filter(|&k| {
                    mutations[*k as usize].position() >= lastpos
                        && mutations[*k as usize].position() < *pos
                })
                .for_each(|a| output.push(*a));
            lastpos = *pos;
            std::mem::swap(&mut parent1_genome, &mut parent2_genome);
        }
        parent1_genome
            .mutations
            .iter()
            .filter(|&k| mutations[*k as usize].position() >= lastpos)
            .for_each(|a| output.push(*a));
        for m in new_mutations {
            output.push(*m);
        }
        output.sort_by(|i, j| {
            mutations[*i as usize]
                .position()
                .cmp(&mutations[*j as usize].position())
        });
        output
    }

    fn keys_to_positions(mutations: &[Mutation], keys: &[u32]) -> Vec<i64> {
        let mut rv = vec![];

        for k in keys {
            rv.push(mutations[*k as usize].position().into())
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
            .all(|w| mutations[w[0] as usize].position() <= mutations[w[1] as usize].position()),);
        assert!(parent2_genome
            .mutations
            .windows(2)
            .all(|w| mutations[w[0] as usize].position() <= mutations[w[1] as usize].position()),);
        let mut offspring_genomes = Vec::<u32>::new();
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
            .all(|w| mutations[w[0] as usize].position() <= mutations[w[1] as usize].position()),);
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
