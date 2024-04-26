To run studies and process data, simply run the `main.m` script. This requires the COMSOL Wave Optics module.    
To just process the provided data without running any simulations, comment out the `Run COMSOL` portion of the `main.m` script and 
- for the study of plane waves normally incident from below the funnel set
  - `study.sim = 0`
  - `study.sweepName = 'results/HMM-sweep'`
- for the study of in-plane polarized dipole radiation at the funnel tip set
  - `study.sim = 1`
  - `study.sweepName = 'results/HMM-sweep-X'`    

then run `main.m`.    
Alternatively, you can call the `processStudy(path)` or `processEmission(path)` functions on incident waves or emitted radiation studies respectively, where `path` is the path to the directory containing the generated data. 
