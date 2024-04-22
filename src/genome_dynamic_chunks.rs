use forrustts::genetics::GenerateBreakpoints;
use forrustts::genetics::GeneticMap;
use forrustts::prelude::*;

use crate::common::DiploidGenome;
use crate::common::Mutation;

struct Chunk {
    mutations: Vec<usize>,
    num_mutations: usize,
    count: u32,
}

#[derive(Default, Clone)]
struct Genome {
    chunks: Vec<usize>,
}

pub struct DiploidPopulation {
    chunk_length: usize,
    chunks: Vec<Chunk>,
    mutations: Vec<Mutation>,
    mutation_counts: Vec<u32>,
    individuals: Vec<DiploidGenome>,
    genomes: Vec<Genome>,
}

impl DiploidPopulation {
    pub fn new(size: usize, chunk_length: usize) -> Self {
        let genomes = vec![Genome::default(); 2 * size];
        let individuals = vec![
            DiploidGenome {
                first: 0,
                second: 0
            };
            size
        ];

        Self {
            chunk_length,
            genomes,
            individuals,
            mutations: vec![],
            mutation_counts: vec![],
            chunks: vec![],
        }
    }
}

fn crossover(
    genome1: usize,
    genome2: usize,
    mutations: &Vec<Mutation>,
    crossover_position: Vec<Position>,
    new_mutations: Vec<usize>,
    genomes: &mut Vec<Genome>,
    chunks: &mut Vec<Chunk>,
) -> Genome {
    let mut current_genome = genome1;
    let mut current_genome_index = 0_usize;
    let mut current_genome_index_position = 0_usize;

    let mut new_chunks = Vec::<usize>::default();

    todo!("not done");

    Genome { chunks: new_chunks }
}
