def conf_generator(smile, numConf, params=None, getParams=None):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdDistGeom

    if getParams is None:
        getParams = False

    # Set ETKDG parameters
    if params is None:
        #params = rdDistGeom.srETKDGv3()
        params = rdDistGeom.ETKDGv3()
        params.useRandomCoords = True
        params.numThreads = 0
        params.maxAttempts = 10000
        params.pruneRmsThresh = 1.0

    if getParams:
        parameters = dir(params)
        parameters = [x for x in parameters if '__' not in x]
        parameters = [x for x in parameters if 'Set' not in x]
        parameter_val = [getattr(params, x) for x in parameters]
        param_dict = dict(zip(parameters, parameter_val))
        return param_dict

    # Create molecule
    mol = Chem.MolFromSmiles(smile)
    Chem.SanitizeMol(mol)

    # Add H's
    mol = Chem.AddHs(mol, addCoords=True, explicitOnly=True) #, addResidueInfo=True)


    # Produce conformers
    AllChem.EmbedMultipleConfs(mol, numConfs=numConf, params=params)

    return mol,params
