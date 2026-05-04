def launch_stdep(root_dir: str = './',
                ucell=None,
                make_supercell: bool = False,
                scell_mat=None,
                scell=None,
                temperature: float = None,
                preexisting_ifcs: bool = False,
                preexisting_ifcs_path=None,
                max_freq: float = None,
                quantum: bool = True,
                tdep_bin_directory: str = None,
                first_order: bool = True,
                displ_threshold_firstorder: float = 0.0001,
                max_iterations_first_order: int = None,
                rc2: float = None,
                polar: bool = False,
                loto_infile: str = None,
                niters: int = 10,
                nconfs: List[int] = None,
                remove_infiles: bool = True,
                restart = False,
                mlip_bin: str = './mlp',
                bin_prefix: str = '',
                mlip_model=None): ##### !! MlipModel object !! ####
    
    # sTDEP steps:
    # 1. first iteration
    #     1.1 sample 4 configurations with the initial ifcs
    #     1.2 compute properties 
    # 2. loop over iterations
    #     2.1 extract ifcs (with first order) from previous trajectories
    #     2.2 sample N(i) configurations with the previous ifcs
    #     2.3 compute properties
    # 3. final step
    #    3.1 extract ifcs with the final iteration configurations

    if make_supercell == False:
        if scell is None:
            raise TypeError('Since make_supercell is False, you must provide a supercell')
    else:
        if scell_mat is None:
            raise TypeError('Since make_supercell is True, you must provide scell_mat!')
        else:
            scell = mk_supercell(ucell, scell_mat) 
    if preexisting_ifcs is False:
        if max_freq is None:
            raise TypeError('Since preexisting_ifcs is False, you must provide max_freq!')
    else:
        preexisting_ifcs_path = Path(preexisting_ifcs_path)
        if not preexisting_ifcs_path.is_file():
            raise ValueError(f'The file {preexisting_ifcs_path.absolute()} does not exist!')
        max_freq = False
    if tdep_bin_directory is not None:
        tdep_bin_directory = Path(tdep_bin_directory)
    
    if polar == True:
        if loto_infile is None:
            raise TypeError('Since polar is True, you must provide loto_infile!')
        else:
            loto_infile = Path(loto_infile)
    
    if niters != len(nconfs):
        raise ValueError(f'The length of `nconfs` must be equal to the number of iterations (`niters`)')
    niters = len(nconfs)

    root_dir = Path(root_dir)

    # infiles
    infiles_dir = root_dir.joinpath('infiles')
    infiles_dir.mkdir(parents=True, exist_ok=True)
    write(infiles_dir.joinpath('infile.ucposcar'), ucell, format='vasp')
    ucell_path = infiles_dir.joinpath('infile.ucposcar')
    scell_path = infiles_dir.joinpath('infile.ssposcar')
    write(infiles_dir.joinpath('infile.ssposcar'), scell, format='vasp')

    if preexisting_ifcs == True:
        shutil.copy(preexisting_ifcs_path, infiles_dir.joinpath('infile.forceconstant'))

    iters_dir = root_dir.joinpath('iterations')
    iters_dir.mkdir(parents=True, exist_ok=True)

    print_kb('++++++++++++++++++++++++++++++++++++++')
    print_kb('----------- sTDEP launched -----------')
    print_kb('++++++++++++++++++++++++++++++++++++++')
    results = []

    # check if it's a restart
    if restart == True:
        prev_iter_dirs = [x for x in iters_dir.glob('iter_*')]
        prev_iters = [int(x.name.split('_')[-1]) for x in prev_iter_dirs]
        done_prev_iters = []
        undone_prev_iters = []
        for prev_iter in prev_iters:
            iter_dir = iters_dir.joinpath(f'iter_{prev_iter}')
            confs_dir = iter_dir.joinpath('true_props/new_confs_computed.traj')
            if confs_dir.is_file():
                done_prev_iters.append(prev_iter)
            else:
                undone_prev_iters.append(prev_iter)
        del prev_iters
        if len(done_prev_iters) == 0:
            raise ValueError('You asked to restart, but no previous complete iteration was found!')
        else:
            last_iter_done = max(done_prev_iters)
            undone_prev_iters_to_delete = [x for x in undone_prev_iters if x > last_iter_done] # this might be redundant, since undone_prev_iters_to_delete == undone_prev_iters (in principles) but I keep it
            print(f'You asked to restart: the last completed iteration is n. {last_iter_done}; any successive non-complete iteration will be deleted.')
            [shutil.rmtree(iters_dir.joinpath(f'iter_{x}')) for x in undone_prev_iters_to_delete]
            first_iteration = last_iter_done + 1
    else:
        first_iteration = 1
        
    # 1. first iteration
    for iter in range(first_iteration, niters+1):
        print(f'====== ITERATION n. {iter} ======')
        iter_dir = iters_dir.joinpath(f'iter_{iter}')
        confs_dir = iter_dir.joinpath('configurations')
        confs_dir.mkdir(parents=True, exist_ok=True)

        #   1.1 sample 4 configurations with the initial ifcs

        if iter != 1:
            ifc_dir = iter_dir.joinpath('ifcs')
            ifc_dir.mkdir(parents=True, exist_ok=True)
            print('I after first')
            # extract force constants from previous confs
            previous_iter_dir = iters_dir.joinpath(f'iter_{iter-1}')
            # we need the previous confs
            traj_to_get_path = previous_iter_dir.joinpath('true_props/new_confs_computed.traj')
            #print(previous_iter_dir.absolute(), traj_to_get_path.absolute())
            traj = read(traj_to_get_path.absolute(), index=':')
            # and the previous unitcell + supercell
            ucell_path = previous_iter_dir.joinpath('true_props/infile.ucposcar')
            scell_path = previous_iter_dir.joinpath('true_props/infile.ssposcar')
            ucell = read(ucell_path, format='vasp')
            scell = read(scell_path, format='vasp')

            # extract ifcs (only 2nd order)
            tdp.extract_ifcs(from_infiles = False,
                            infiles_dir = None,
                            unitcell = ucell,
                            supercell = scell,
                            sampling = traj,
                            timestep = 1,
                            dir = ifc_dir.absolute(),
                            first_order = first_order,
                            displ_threshold_firstorder = displ_threshold_firstorder,
                            max_iterations_first_order = max_iterations_first_order,
                            rc2 = rc2, 
                            rc3 = None, 
                            polar = polar,
                            loto_filepath = loto_infile, 
                            stride = 1, 
                            temperature = temperature,
                            bin_prefix = bin_prefix,
                            tdep_bin_directory = tdep_bin_directory)
            
            if first_order == True:
                ucell_path = ifc_dir.joinpath('first_order_optimisation/optimized_unitcell.poscar')
                scell_path = ifc_dir.joinpath('first_order_optimisation/optimized_supercell.poscar')
                ucell = read(ucell_path)
                scell = read(scell_path)
            
        
            last_ifc_path = ifc_dir.joinpath('outfile.forceconstant')

            if remove_infiles == True:
                names = ['stat', 'meta', 'positions', 'forces']
                [ifc_dir.joinpath(f'infile.{name}').unlink() for name in names]
                if first_order == True:
                    [ifc_dir.joinpath(f'first_order_optimisation/infile.{name}').unlink() for name in names]
        else:
            # it's the first iteration
            pass

        # now ucell is the unitcell used for the extraction (after first order, if done), same for the scell 
        make_canonical_configurations_parameters = dict(ucell = ucell,
                                                        scell = scell,
                                                        nconf = nconfs[iter-1],
                                                        temp = temperature,
                                                        quantum = quantum,
                                                        dir = confs_dir,
                                                        outfile_name = 'new_confs.traj', # this will be saved inside dir
                                                        pref_bin = bin_prefix,
                                                        tdep_bin_directory=tdep_bin_directory)
             
        if iter == 1:
            if preexisting_ifcs == False:
                make_canonical_configurations_parameters['max_freq'] = max_freq
            else:
                make_canonical_configurations_parameters['ifcfile_path'] = preexisting_ifcs_path
        else:
            make_canonical_configurations_parameters['ifcfile_path'] = last_ifc_path 
        
        tdp.make_canonical_configurations(**make_canonical_configurations_parameters)

        latest_confs = read(confs_dir.joinpath('new_confs.traj'), index=':')

        # compute properties
        prop_iter = iter_dir.joinpath('true_props')
        prop_iter.mkdir(exist_ok=True, parents=True)

        parameters = {'mlip_bin' : mlip_bin, 'bin_pref' : bin_prefix}
        new_confs_computed = mlip_model.compute_properties(atoms=latest_confs,
                                                           wdir=prop_iter,
                                                           parameters=parameters)
        write(prop_iter.joinpath('new_confs_computed.traj'), new_confs_computed)
       
        ln_s_f(ucell_path, prop_iter.joinpath('infile.ucposcar'))
        ln_s_f(scell_path, prop_iter.joinpath('infile.ssposcar'))
    

    