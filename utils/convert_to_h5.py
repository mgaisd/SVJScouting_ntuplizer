from __future__ import print_function, division
import os
import argparse
import uproot
import numpy as np
import h5py
import awkward as ak

def to_np_array(ak_array, maxN=100, pad=0):
    '''convert awkward array to regular numpy array'''
    return ak.to_numpy(ak.fill_none(ak.pad_none(ak_array,maxN,clip=True,axis=-1),pad))

def store_objects_coordinates(arrays, nentries, nobj=10, obj='FatJet_'):
    '''store objects in zero-padded numpy arrays'''
    l1Obj_cyl = np.zeros((nentries,nobj,3))
    l1Obj_cart = np.zeros((nentries,nobj,3))
    pt = to_np_array(arrays['{}pt'.format(obj)],maxN=nobj)
    eta = to_np_array(arrays['{}eta'.format(obj)],maxN=nobj)
    phi = to_np_array(arrays['{}phi'.format(obj)],maxN=nobj)
    l1Obj_cyl[:,:,0] = pt
    l1Obj_cyl[:,:,1] = eta
    l1Obj_cyl[:,:,2] = phi
    l1Obj_cart[:,:,0] = pt*np.cos(phi)
    l1Obj_cart[:,:,1] = pt*np.sin(phi)
    l1Obj_cart[:,:,2] = pt*np.sinh(eta)
    
    return l1Obj_cyl, l1Obj_cart

def store_objects_truth(arrays, nentries, nobj=10, obj='FatJet_'):
    '''store objects in zero-padded numpy arrays'''
    l1Obj_truth = np.zeros((nentries,nobj,1))
    fromsuep = to_np_array(arrays['{}fromsuep'.format(obj)],maxN=nobj)
    l1Obj_truth[:,:,0] = fromsuep
    
    return l1Obj_truth

def store_objects_addfeatures(arrays, nentries, nobj=10, obj='FatJet_'):
    '''store objects in zero-padded numpy arrays'''
    l1Obj_features = np.zeros((nentries,nobj,2))
    pdgid = to_np_array(arrays['{}pdgid'.format(obj)],maxN=nobj)
    fjidx = to_np_array(arrays['{}fjidx'.format(obj)],maxN=nobj)
    l1Obj_features[:,:,0] = pdgid
    l1Obj_features[:,:,1] = fjidx
    
    return l1Obj_features

def store_objects_features(arrays, nentries, nobj=10,obj='FatJet_'):
    '''store objects in zero-padded numpy arrays'''
    features=arrays.fields
    nfeats = len(arrays.fields)
    objects = np.zeros((nentries,nobj,nfeats+3)) #+3 for cart. coordinates
    for i in range(0,nfeats):
        objects[:,:,i] = to_np_array(arrays['{}'.format(features[i])],maxN=nobj)
    #add cartesian coordinates
    pt = to_np_array(arrays['{}_pt'.format(obj)],maxN=nobj)
    eta = to_np_array(arrays['{}_eta'.format(obj)],maxN=nobj)
    phi = to_np_array(arrays['{}_phi'.format(obj)],maxN=nobj)
    objects[:,:,nfeats] = pt*np.cos(phi) #px
    objects[:,:,nfeats+1] = pt*np.sin(phi) #py
    objects[:,:,nfeats+2] = pt*np.sinh(eta) #pz
    features_names = features + ['px','py','pz']
    features_names = [n.encode('utf8') for n in features_names]
    return objects, features_names

def convert_event_based(input_file, output_file, tree_name):
    inFile = uproot.open(input_file)
    l1Tree = inFile[tree_name]
    nentries = l1Tree.num_entries

    # save up to n objects (jets, muons, electrons)
    njets = 20
    nfatjets = 10
    nmuons = 6
    nelectrons = 6
    nphotons = 20
    npfcands=1000

    
    cylNames = [b'pt', b'eta', b'phi']
    cartNames = [b'px', b'py', b'pz']

    # variables to retrieve
    varList = ['n_jet','n_fatjet','n_pho',
               'n_ele', 'n_mu', 'n_pfcand', 'n_bpfcand']
    common_prop = ['pt','eta','phi']
    varList += ['Electron_'+p for p in common_prop+['m']]
    varList += ['Photon_'+p for p in common_prop+['m']]
    varList += ['Muon_'+p for p in common_prop+['m']]
    varList += ['Jet_'+p for p in common_prop+['m']]
    varList += ['FatJet_'+p for p in common_prop+['mass','msoftdrop', 'mtrim']]
    varList += ['PFcand_'+p for p in common_prop+['m','pdgid','fjidx','fromsuep']]
    varList += ['bPFcand_'+p for p in common_prop+['m']]


    # get awkward arrays
    arrays = l1Tree.arrays(varList)

    # store objects: jets, muons, electrons
    l1Jet_cyl, l1Jet_cart = store_objects_coordinates(arrays, nentries, nobj=njets, obj='Jet_')
    l1FatJet_cyl, l1FatJet_cart = store_objects_coordinates(arrays, nentries, nobj=nfatjets, obj='FatJet_')
    l1mu_cyl, l1mu_cart = store_objects_coordinates(arrays, nentries, nobj=nmuons, obj='Muon_')
    l1pho_cyl, l1pho_cart = store_objects_coordinates(arrays, nentries, nobj=nphotons, obj='Photon_')
    l1ele_cyl, l1ele_cart = store_objects_coordinates(arrays, nentries, nobj=nelectrons, obj='Electron_')
    l1pfcand_cyl, l1pfcand_cart = store_objects_coordinates(arrays, nentries, nobj=npfcands, obj='PFcand_')
    l1pfcand_feat = store_objects_addfeatures(arrays, nentries, nobj=npfcands, obj='PFcand_')
    l1pfcand_truth = store_objects_truth(arrays, nentries, nobj=npfcands, obj='PFcand_')
    l1bpfcand_cyl, l1bpfcand_cart = store_objects_coordinates(arrays, nentries, nobj=npfcands, obj='bPFcand_')

    
    with h5py.File(output_file, 'w') as outFile:
        outFile.create_dataset('FeatureNames_cyl', data=cylNames, compression='gzip')
        outFile.create_dataset('FeatureNames_cart', data=cartNames, compression='gzip')
        outFile.create_dataset('Jet_cyl', data=l1Jet_cyl, compression='gzip')
        outFile.create_dataset('Jet_cart', data=l1Jet_cart, compression='gzip')
        outFile.create_dataset('FatJet_cyl', data=l1FatJet_cyl, compression='gzip')
        outFile.create_dataset('FatJet_cart', data=l1FatJet_cart, compression='gzip')
        outFile.create_dataset('Muon_cyl', data=l1mu_cyl, compression='gzip')
        outFile.create_dataset('Muon_cart', data=l1mu_cart, compression='gzip')
        outFile.create_dataset('Ele_cyl', data=l1ele_cyl, compression='gzip')
        outFile.create_dataset('Ele_cart', data=l1ele_cart, compression='gzip')
        outFile.create_dataset('Pho_cyl', data=l1pho_cyl, compression='gzip')
        outFile.create_dataset('Pho_cart', data=l1pho_cart, compression='gzip')
        outFile.create_dataset('Pfcand_cyl', data=l1pfcand_cyl, compression='gzip')
        outFile.create_dataset('Pfcand_cart', data=l1pfcand_cart, compression='gzip')
        outFile.create_dataset('Pfcand_feat', data=l1pfcand_feat, compression='gzip')
        outFile.create_dataset('Pfcand_truth', data=l1pfcand_truth, compression='gzip')
        outFile.create_dataset('bPfcand_cyl', data=l1bpfcand_cyl, compression='gzip')
        outFile.create_dataset('bPfcand_cart', data=l1bpfcand_cart, compression='gzip')

def convert_jet_based(input_file, output_file, tree_name):
    inFile = uproot.open(input_file)
    l1Tree = inFile[tree_name]
    nentries = l1Tree.num_entries

    jetFeatureNames = ['area','n2b1','n3b1','tau1','tau2','tau3','tau4','mass','msoftdrop','mtrim']
    particleFeatureNames = ['pdgid','fjidx'] 
    common_prop = ['pt','eta','phi']
    jetFeatureNames += [p for p in common_prop]
    particleFeatureNames += [p for p in common_prop]

    # variables to retrieve
    varList = ['n_pfcand','n_fatjet']
    varJets = ['FatJet_'+p for p in jetFeatureNames]
    varJets += ['n_fatjet']
    varPfcands = ['PFcand_'+p for p in particleFeatureNames+['m']]
    varList+=varJets
    varList+=varPfcands
    
    # get awkward arrays
    arrays = l1Tree.arrays(varList)    

    jet_record = ak.zip({"FatJets": ak.zip({ name : arrays[name] for name in varJets})})#
    pfcand_record = ak.zip({"pfcands": ak.zip({ name : arrays[name] for name in varPfcands})})
    # save up to 1 fat jet for simplicity
    nfatjets = 1
    #remove empty fat jet arrays 
    #jet_record.FatJets = jet_record.FatJets[jet_record.FatJets.n_fatjet > 0] # or we can do pt>0 , need to check the dataset
    #jet_record.FatJets = jet_record.FatJets[jet_record.FatJets.FatJet_pt > 0] # or we can do pt>0 , need to check the dataset
    #choose pf candidates associated to this one fat jet
    fatjetidx = 0
    pfcand_record.pfcands = pfcand_record.pfcands[pfcand_record.pfcands.PFcand_fjidx == fatjetidx]
    
    # store objects: jets, and pfcands
    fatjets,fatjets_names = store_objects_features(jet_record.FatJets, nentries, nobj=nfatjets,obj='FatJet')
    pfcands,pfcands_names = store_objects_features(pfcand_record.pfcands, nentries, nobj=100,obj='PFcand')

    
    with h5py.File(output_file, 'w') as outFile:
        outFile.create_dataset('jetFeatureNames', data=fatjets_names, compression='gzip')
        outFile.create_dataset('particleFeatureNames', data=pfcands_names, compression='gzip')
        outFile.create_dataset('fatjets', data=fatjets, compression='gzip')
        outFile.create_dataset('jetConstituentList', data=pfcands, compression='gzip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outtype', type=str, required=True,help='Use event or jet for event/jet based output')
    parser.add_argument('--inpfile', type=str, required=True)
    parser.add_argument('--outfile', type=str, required=True)
    parser.add_argument('--treename', type=str, default='mmtree/tree')
    args = parser.parse_args()
    if args.outtype=='event' :  
        convert_event_based(args.inpfile, args.outfile, args.treename)
    elif args.outtype=='jet' :
        convert_jet_based(args.inpfile, args.outfile, args.treename)

