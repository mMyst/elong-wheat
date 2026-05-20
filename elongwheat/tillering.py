from openalea.plantgl.algo import *
import math


from elongwheat import parameters as elongwheat_parameters

def update_tiller_replications(g, adel_wheat, plant_density,tillers_replications, gaic, coef_buffer_til, coef_delay_til,
                               GAIp=None):
    """
    Evaluates thermal time since primordium and updates tiller replications 
    if the conditions for tiller emergence are met.
    """
    # Fetch plastochrone from parameters
    plastochrone = elongwheat_parameters.PARAMETERS.PLASTOCHRONE
    ms_vid = next(vid for vid, label in g.property('label').items() if label == 'MS') #get MS vid id 

    for vid in g.components_at_scale(ms_vid,3):
        tiller_rank = g.index(vid)   
        if (g.label(vid) == 'metamer0' or 'T'+tiller_rank in tillers_replications.keys()) or 'hiddenzone' not in g.get_vertex_property(vid).keys() :
            continue

        hz_age = g.get_vertex_property(vid)['hiddenzone']['hiddenzone_age']
        

        in_window = (coef_delay_til*plastochrone <= hz_age < plastochrone*(coef_delay_til +coef_buffer_til))
        
        if in_window :
            if GAIp is None:
                GAIp = compute_1plant_GAIp(adel_wheat.scene(g), plant_density)   
            
            print(f" 🚨 Simulation gaic_{gaic}_dens_{plant_density}_delay_{coef_delay_til}_buf_{coef_buffer_til} is in window")    

            tillers_replications = tiller_initiation(
                tillers_replications, 
                tiller_rank, 
                GAIp, 
                GAIc=gaic
            ) 
            
    return tillers_replications


def tiller_initiation(tillers_replications, tiller_rank, GAIp,GAIc=0.6):
    """
    decide tiller emission based on GAI.
    :param ??? initiated_leaves: number of initiated leaves on the MS, used to determine the rank of the initiating tiller 
    rq : puisque plastochrone est fixe (?) on pourrait appeler la fonction seulement autour des TpsThq modulo plastochrone = 0 
    c'est ce qu'on disait de faire avec teq_since_primordium <= conditionne l'appel à cette fonction

    :param float GAIp: current GAI of the plant
    :param float GAIc: GAI threshold for tiller emergence. default is 0.6, from Lecarpentier et al. 2019 (WALTer)
    :param float delta_til: time window for tiller emergence 
        (thermal time. 500DD by default. should otherwise be calculated as a fraction of the phyllochron) 
    :param dict tiller_replication: dictionnary with metabolic wheight of each tiller. starts as none or {} when there is no tillers
    :return: tiller_replication e.g. : {'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}
    """
    
    print('GAIp is ' + str(GAIp))
    print("Tiller #" + str(tiller_rank) + " is trying to initiate...")
    

    #if  -delta_til < teq_since_primordium < delta_til : #buffer window during which initiation is possible => move up
    if GAIp < GAIc : #on/off switch. could be reworked into a probability function
        # if tillers_replications is None: #a voir si on décide d'initialiser tiller replication comme None ou comme {}
        #     tillers_replications = {}
        tillers_replications["T"+tiller_rank] = 0.5
        print("...and it did !")
    else:
        
        print("...but it failed :(")


    return tillers_replications


def compute_1plant_GAIp(scene_1plante,density):
    """Compute the GAI of a one plante scene
    Parameters
    ----------
    scene_1plante : scene of the plant
    density : plant density (plant/m2)
    """

    GA = openalea.plantgl.algo.surface(scene_1plante)
    GAIp = GA * density
    return GAIp

def compute_scene_GAIp(scene,domain):
    """Compute the GAI of a 3D scene
    Parameters
    ----------
    scene : scene of the plant
    """

    GA = openalea.plantgl.algo.surface(scene)
    GAIp = GA / domain
    return GAIp
