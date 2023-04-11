use rand::prelude::Rng;
use rand::prelude::SeedableRng;

use crate::common::generate_mutations;
use crate::common::DiploidGenome;
use crate::common::Mutation;
use crate::common::ParentalGenome;
use crate::SimParams;

use forrustts::genetics::GenerateBreakpoints;
use forrustts::genetics::GeneticMap;
use forrustts::prelude::*;

#[derive(Default)]
struct HaploidGenome {
    mutations: Vec<usize>,
    count: u32,
}

// The borrow checker is not gonna like this...
pub struct DiploidPopulation {
    haplotypes: Vec<HaploidGenome>,
    individuals: Vec<DiploidGenome>,
    mutations: Vec<Mutation>,
    mutation_counts: Vec<u32>,
}

impl DiploidPopulation {
    pub fn new(size: usize) -> Self {
        let haplotypes = vec![HaploidGenome {
            mutations: vec![],
            count: 2 * (size as u32),
        }];
        let individuals = vec![
            DiploidGenome {
                first: 0,
                second: 0
            };
            size
        ];

        Self {
            haplotypes,
            individuals,
            mutations: vec![],
            mutation_counts: vec![],
        }
    }

    pub fn num_segregating_mutations(&self) -> u32 {
        todo!("oops, still need this");
    }
}

fn make_haploid_genome_queue(genomes: &[HaploidGenome]) -> Vec<usize> {
    let mut rv = vec![];
    genomes.iter().enumerate().for_each(|(i, g)| {
        if g.count == 0 {
            rv.push(i);
        }
    });
    rv
}

fn get_parental_genome(genomes: &[HaploidGenome], genome: usize) -> ParentalGenome {
    ParentalGenome {
        mutations: &genomes[genome].mutations,
        current_mutation_index: 0,
        genome,
    }
}

fn get_parental_genomes<'a>(
    parents: &'a [DiploidGenome],
    genomes: &'a [HaploidGenome],
    parent: usize,
) -> (ParentalGenome<'a>, ParentalGenome<'a>) {
    (
        get_parental_genome(genomes, parents[parent].first),
        get_parental_genome(genomes, parents[parent].second),
    )
}

#[inline(never)]
pub fn evolve_pop_with_haplotypes(
    params: SimParams,
    genetic_map: GeneticMap,
) -> Option<DiploidPopulation> {
    let params = params.validate()?;
    let mut pop = DiploidPopulation::new(params.num_individuals as usize);

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
    for generation in 0..params.num_generations {
        let mut offspring: Vec<DiploidGenome> = vec![];
        let mut genome_queue = make_haploid_genome_queue(&pop.haplotypes);
        for _ in 0..params.num_individuals {
            // Pick two parents
            let parent1 = rng.sample(parent_picker);
            let parent2 = rng.sample(parent_picker);

            let genomes = get_parental_genomes(&pop.individuals, &pop.haplotypes, parent1);
            let genomes = if rng.sample(u01) < 0.5 {
                genomes
            } else {
                (genomes.1, genomes.0)
            };

            // Mutations for offspring genome 1
            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut vec![],
                &mut pop.mutations,
                &mut rng,
            );

            genetic_map.generate_breakpoints(&mut rng);
            if mutations.is_empty() && genetic_map.breakpoints().is_empty() {
                pop.haplotypes[genomes.0.genome].count += 1;
            } else {
                unimplemented!("get a new genome and call our lib fn to populate it");
            }
        }
    }
    Some(pop)
}
