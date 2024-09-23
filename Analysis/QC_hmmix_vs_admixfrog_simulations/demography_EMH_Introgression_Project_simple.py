from collections import defaultdict, namedtuple

import numpy as np
import msprime
import pandas as pd
import stdpopsim



#helper functions
MERGE = namedtuple("merge", ('source', 'target', 'time'))



# msprime version > 1
def define_popconfig(pop_params):
    """Generate list of population configurations for msprime."""
    demographic_events = msprime.Demography()
    for p in pop_params:
        demographic_events.add_population(name=str(p),initial_size=pop_params[p]["Ne"], initially_active=True)
    return demographic_events


def define_samples(samples,generation_time):
    """Generate list of sample definitions for msprime."""
    g_time = generation_time/1000
    sample_names = []
    for i, pop in enumerate(samples):
        times = [int(t/g_time) for t in samples[pop]["t_sample"]]
        sample_names.extend([msprime.SampleSet(1,population=i, time=int(t),ploidy=1) for t in times])
    return sample_names


def demo_EMH_Introgression_Project_simplest(
    generation_time,
    admixture_proportion,
    print_Demo=False
    ):

    generation_time = generation_time/1000

    #population sizes, size changes and labels
    n_N = 1000
    n_D = 1000
    n_Hafr = 10000
    n_Heur = 10000
    n_Hasn = 10000
    n_P = 10000


    tm_NHeur = 50
    m_NHeur = float(admixture_proportion)


    pop_params = {
        "nea":  {"id": 0, "Ne": n_N},
        "den":  {"id": 1, "Ne": n_D},
        "afr": {"id": 2, "Ne": n_Hafr },
        #"eur": {"id": 3, "Ne": n_Heur, "t_bottle": t_exp, "Ne_bottle": n_ancHeur,"t_OoA_bottle": (70 - 0.1)/generation_time, "Ne_OoA_bottle": 250},
        #"asn": {"id": 4, "Ne": n_Hasn, "t_bottle": t_exp, "Ne_bottle": n_ancHeur},
        "eur": {"id": 3, "Ne": n_Heur},
        "asn": {"id": 4, "Ne": n_Hasn},
        "cont": {"id": 5, "Ne": n_Heur},
        "pan": {"id": 6, "Ne": n_P}
    }
    ids = dict((k, p["id"]) for k, p in pop_params.items())

    demographic_events = define_popconfig(pop_params)


    """3. set up base tree"""
    base_tree = [  #give times in years
        MERGE('cont' ,'eur', 0+(generation_time)*2),
        MERGE('asn' ,'eur', 45),
        MERGE('eur' ,'afr', 70),
        MERGE('nea' ,'den', 400),
        MERGE('afr' ,'den', 600),
        MERGE('pan' ,'den', 6000),
    ]


    demographic_events.add_migration_rate_change(time=int(tm_NHeur/generation_time),rate= m_NHeur,source=3,dest=0),
    demographic_events.add_migration_rate_change(time=int(tm_NHeur/generation_time)+1,rate=0,source=3,dest=0)

    for m in base_tree:
        demographic_events.add_population_split(time=(m.time / generation_time), derived=[str(m.source)], ancestral=str(m.target))

    demographic_events.sort_events()

    if print_Demo == True:
        print("Printing demography")
        check=demographic_events.debug()
        check.print_history(output=open("../EMH_Introgression_Project_sim_demo_simplest.txt", 'w'))

    return demographic_events


def demo_EMH_Introgression_Project_simplified(
    generation_time,
    admixture_proportion,
    print_Demo=False
    ):

    generation_time = generation_time/1000

    #population sizes, size changes and labels
    n_N = 1000
    n_D = 1000
    n_Hafr = 10000
    n_Heur = 10000
    n_Hasn = 10000
    n_P = 10000

    t_exp = 10/generation_time
    n_ancHeur = 2000

    tm_NHeur = 50
    m_NHeur = float(admixture_proportion)


    pop_params = {
        "nea":  {"id": 0, "Ne": n_N},
        "den":  {"id": 1, "Ne": n_D},
        "afr": {"id": 2, "Ne": n_Hafr },
        "eur": {"id": 3, "Ne": n_Heur, "t_OoA_bottle": (70 - 0.1)/generation_time, "Ne_OoA_bottle": 250},
        "asn": {"id": 4, "Ne": n_Hasn},
        "cont": {"id": 5, "Ne": n_Heur},
        "pan": {"id": 6, "Ne": n_P}
    }
    ids = dict((k, p["id"]) for k, p in pop_params.items())

    demographic_events = define_popconfig(pop_params)


    # add stron OoA bottleneck from 69900 to 70000 kya
    demographic_events.add_population_parameters_change(
        time= pop_params["eur"]["t_OoA_bottle"],
        initial_size=pop_params["eur"]["Ne_OoA_bottle"],
        population=ids['eur'])




    """3. set up base tree"""
    base_tree = [  #give times in years
        MERGE('cont' ,'eur', 0+(generation_time)*2),
        MERGE('asn' ,'eur', 45),
        MERGE('eur' ,'afr', 70),
        MERGE('nea' ,'den', 400),
        MERGE('afr' ,'den', 600),
        MERGE('pan' ,'den', 6000),
    ]


    demographic_events.add_migration_rate_change(time=int(tm_NHeur/generation_time),rate= m_NHeur,source=3,dest=0),
    demographic_events.add_migration_rate_change(time=int(tm_NHeur/generation_time)+1,rate=0,source=3,dest=0)

    for m in base_tree:
        demographic_events.add_population_split(time=(m.time / generation_time), derived=[str(m.source)], ancestral=str(m.target))

    demographic_events.sort_events()

    if print_Demo == True:
        print("Printing demography")
        check=demographic_events.debug()
        check.print_history(output=open("../EMH_Introgression_Project_simplified_sim_demo.txt", 'w'))

    return demographic_events

def demo_EMH_Introgression_Project_with_diverged_DEN(
    generation_time,
    admixture_proportion,
    admixture_proportion_den,
    print_Demo=False
    ):

    generation_time = generation_time/1000

    #population sizes, size changes and labels
    n_N = 1000
    n_D = 1000
    n_Hafr = 10000
    n_Heur = 10000
    n_Hasn = 10000
    n_P = 10000

    t_exp = 10/generation_time
    n_ancHeur = 2000

    tm_NHeur = 50
    m_NHeur = float(admixture_proportion)

    tm_DHasn = 43
    m_DHasn = float(admixture_proportion_den)


    pop_params = {
        "nea":  {"id": 0, "Ne": n_N},
        "den":  {"id": 1, "Ne": n_D},
        "afr": {"id": 2, "Ne": n_Hafr },
        "eur": {"id": 3, "Ne": n_Heur, "t_OoA_bottle": (70 - 0.1)/generation_time, "Ne_OoA_bottle": 250},
        "asn": {"id": 4, "Ne": n_Hasn},
        "cont": {"id": 5, "Ne": n_Heur},
        "pan": {"id": 6, "Ne": n_P},
        "int_den":  {"id": 7, "Ne": n_D},
    }
    ids = dict((k, p["id"]) for k, p in pop_params.items())

    demographic_events = define_popconfig(pop_params)


    # add stron OoA bottleneck from 69900 to 70000 kya
    demographic_events.add_population_parameters_change(
        time= pop_params["eur"]["t_OoA_bottle"],
        initial_size=pop_params["eur"]["Ne_OoA_bottle"],
        population=ids['eur'])




    """3. set up base tree"""
    base_tree = [  #give times in years
        MERGE('cont' ,'eur', 0+(generation_time)*2),
        MERGE('asn' ,'eur', 45),
        MERGE('eur' ,'afr', 70),
        MERGE('int_den' ,'den', 350),
        MERGE('nea' ,'den', 400),
        MERGE('afr' ,'den', 600),
        MERGE('pan' ,'den', 6000),
    ]


    demographic_events.add_migration_rate_change(time=int(tm_NHeur/generation_time),rate= m_NHeur,source=3,dest=0),
    demographic_events.add_migration_rate_change(time=int(tm_NHeur/generation_time)+1,rate=0,source=3,dest=0)

    demographic_events.add_migration_rate_change(time=int(tm_DHasn/generation_time),rate= m_DHasn,source=3,dest=7),
    demographic_events.add_migration_rate_change(time=int(tm_DHasn/generation_time)+1,rate=0,source=3,dest=7)

    for m in base_tree:
        demographic_events.add_population_split(time=(m.time / generation_time), derived=[str(m.source)], ancestral=str(m.target))

    demographic_events.sort_events()

    if print_Demo == True:
        print("Printing demography")
        check=demographic_events.debug()
        check.print_history(output=open("../EMH_Introgression_Project_simplified_sim_demo.txt", 'w'))

    return demographic_events

