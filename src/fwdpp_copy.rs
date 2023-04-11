use rand::prelude::Rng;
use rand::prelude::SeedableRng;

use crate::common::generate_mutations;
use crate::common::generate_offspring_genome;
use crate::common::DiploidGenome;
use crate::common::Mutation;
use crate::common::ParentalGenome;
use crate::SimParams;

use forrustts::genetics::GenerateBreakpoints;
use forrustts::genetics::GeneticMap;
use forrustts::prelude::*;

#[derive(Default, Debug)]
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

    #[inline(never)]
    fn mutation_recycling(&self) -> Vec<usize> {
        self.mutation_counts
            .iter()
            .enumerate()
            .filter_map(|(index, value)| if value == &0 { Some(index) } else { None })
            .collect::<Vec<usize>>()
    }

    fn count_mutations(&mut self) {
        self.mutation_counts.fill(0);
        self.mutation_counts.resize(self.mutations.len(), 0);
        for g in &self.haplotypes {
            for m in &g.mutations {
                self.mutation_counts[*m] += g.count;
            }
        }
    }

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

fn get_mendelized_parent_genome_indexes(
    individuals: &[DiploidGenome],
    parent: usize,
    u01: rand::distributions::Uniform<f64>,
    rng: &mut rand::rngs::StdRng,
) -> (usize, usize) {
    if rng.sample(u01) < 0.5 {
        (individuals[parent].first, individuals[parent].second)
    } else {
        (individuals[parent].second, individuals[parent].first)
    }
}

#[inline(never)]
fn update_genomes(
    genomes: (usize, usize),
    mutations: Vec<usize>,
    genetic_map: &GeneticMap,
    pop: &mut DiploidPopulation,
    fixation_threshold: u32,
    genome_queue: &mut Vec<usize>,
) -> usize {
    let (genome1, genome2) = genomes;
    let (output_genome_index, update) =
        if mutations.is_empty() && genetic_map.breakpoints().is_empty() {
            pop.haplotypes[genome1].count += 1;
            (genome1, None)
        } else {
            match genome_queue.pop() {
                Some(index) => {
                    pop.haplotypes[index].mutations.clear();
                    assert_eq!(pop.haplotypes[index].count, 0);
                    (
                        index,
                        Some(HaploidGenome {
                            mutations: vec![],
                            count: 0,
                        }),
                    )
                }
                None => {
                    let rv = pop.haplotypes.len();
                    (
                        rv,
                        Some(HaploidGenome {
                            mutations: vec![],
                            count: 0,
                        }),
                    )
                }
            }
        };
    // FIXME: this is super hack-ish:
    // We need to avoid making a fake genome
    // for the case where the genome is being recycled.
    if let Some(mut genome) = update {
        if output_genome_index < pop.haplotypes.len() {
            std::mem::swap(&mut pop.haplotypes[output_genome_index], &mut genome);
            genome.mutations.clear();
        }
        let genomes = (
            get_parental_genome(&pop.haplotypes, genome1),
            get_parental_genome(&pop.haplotypes, genome2),
        );
        generate_offspring_genome(
            genomes,
            &pop.mutations,
            &pop.mutation_counts,
            fixation_threshold,
            mutations,
            genetic_map.breakpoints(),
            &mut genome.mutations,
        );
        genome.count += 1;
        if output_genome_index < pop.haplotypes.len() {
            std::mem::swap(&mut pop.haplotypes[output_genome_index], &mut genome);
        } else {
            pop.haplotypes.push(genome);
        }
    }

    output_genome_index
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
    let mut offspring: Vec<DiploidGenome> = vec![];
    for generation in 0..params.num_generations {
        let mut genome_queue = make_haploid_genome_queue(&pop.haplotypes);
        for g in &mut pop.haplotypes {
            g.count = 0;
        }
        let mut queue = pop.mutation_recycling();
        for _ in 0..params.num_individuals {
            // Pick two parents
            let parent1 = rng.sample(parent_picker);
            let parent2 = rng.sample(parent_picker);

            let (genome1, genome2) =
                get_mendelized_parent_genome_indexes(&pop.individuals, parent1, u01, &mut rng);

            // Mutations for offspring genome 1
            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut queue,
                &mut pop.mutations,
                &mut rng,
            );

            pop.mutation_counts
                .resize(pop.mutation_counts.len() + mutations.len(), 0);

            genetic_map.generate_breakpoints(&mut rng);
            let nm = mutations.len();
            let first = update_genomes(
                (genome1, genome2),
                mutations,
                &genetic_map,
                &mut pop,
                2 * params.num_individuals,
                &mut genome_queue,
            );
            assert!(pop.haplotypes[first].count > 0);
            if nm > 0 {
                assert!(
                    !pop.haplotypes[first].mutations.is_empty(),
                    "{first} {:?}, {:?}",
                    pop.haplotypes[first],
                    nm,
                );
            }

            let (genome1, genome2) =
                get_mendelized_parent_genome_indexes(&pop.individuals, parent2, u01, &mut rng);

            // Mutations for offspring genome 2
            let mutations = generate_mutations(
                generation,
                num_mutations,
                position_generator,
                &mut queue,
                &mut pop.mutations,
                &mut rng,
            );
            pop.mutation_counts
                .resize(pop.mutation_counts.len() + mutations.len(), 0);

            genetic_map.generate_breakpoints(&mut rng);
            let second = update_genomes(
                (genome1, genome2),
                mutations,
                &genetic_map,
                &mut pop,
                2 * params.num_individuals,
                &mut genome_queue,
            );
            assert!(pop.haplotypes[second].count > 0);
            offspring.push(DiploidGenome { first, second });
        }
        std::mem::swap(&mut pop.individuals, &mut offspring);
        //for i in &pop.individuals {
        //    assert!(pop.haplotypes[i.first].count > 0);
        //    assert!(pop.haplotypes[i.second].count > 0);
        //}
        pop.count_mutations();
        offspring.clear();
    }
    Some(pop)
}
