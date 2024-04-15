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
    chunks: Vec<Chunk>,
    mutations: Vec<Mutation>,
    mutation_counts: Vec<u32>,
    individuals: Vec<DiploidGenome>,
    genomes: Vec<Genome>,
}

impl DiploidPopulation {
    pub fn new(size: usize) -> Self {
        let genomes = vec![Genome::default(); 2 * size];
        let individuals = vec![
            DiploidGenome {
                first: 0,
                second: 0
            };
            size
        ];

        Self {
            genomes,
            individuals,
            mutations: vec![],
            mutation_counts: vec![],
            chunks: vec![],
        }
    }
}
