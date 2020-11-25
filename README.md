# robust-lora

## Environments

Python 3.7

[ns3-3.31](https://www.nsnam.org/releases/ns-3-31/) + [lorawan module](https://github.com/signetlabdei/lorawan)

MATLAB 2020a + [SNOPT 7.7](https://ccom.ucsd.edu/~optimizers/solvers/snopt/)

## File Structure

```
.
├── LICENSE
├── README.md    // This file
├── alg          // Algorithms for LoRa gateway placement and device configuration
├── data         // End device location and how to generate path loss matrix 
├── ns3-exp      // Scripts to test in ns3
└── relaxOpt     // MATLAB scripts to call SNOPT to optimally solve the relaxed problem
```

## Algorithms

In `./alg` folder, we implement the following algorithms:

* Proposed greedy gateway placement and device configuration heuristic (`./alg/RGreedy.py`).
* The above proposed heuristic with clustering acceleration (`./alg/clustering.py`).
* Genetic algorithm with the [geneticalgorithm](https://pypi.org/project/geneticalgorithm/) package v1.0.1 (`./alg/RGenetic.py`).

The following baselines are included:

* Greedy heuristic in [Ousat 2019](https://ieeexplore.ieee.org/abstract/document/8815670) (`./alg/ICIOT.py`)

`./alg/main.py` sets which algorithm to run and the parameters of the problem.

To run the algorithms

```bash
python3 ./alg/main.py
```

## ns-3 Experiments

To run the ns-3 simulations, first install [ns3-3.31](https://www.nsnam.org/releases/ns-3-31/).

Then, clone our modified lorawan module and copy the test script

```bash
cd root-of-ns3/ns-3.31/src
git clone https://github.com/Orienfish/lorawan.git
cp path-to-this-repo/ns3-exp/adr.cc root-of-ns3/ns-3.31/scratch
```

To run the ns-3 simulation

```bash
cd root-of-ns3/ns-3.31
./waf --run adr
```

Multiple parameters can be set with the command:

```bash
./waf --run "adr --MType=Confirmed --intfrPowerdBm=-126 --nPeriods=72"
```

For more details, check the help function and the source code.

## Relax Optimization

We develop a relaxed problem and solve it optimally with SNOPT. The implementation is done in MATLAB using interfaces to SNOPT (in `./relaxOpt` folder).

To compare the relaxed bound with the greedy solution under randomly placed end devices, first run algorithms with `python3 ./alg/main.py`, then the generated device locations, candidate gateway locations and necessary communication matrices are directly saved to `./relaxOpt`.

Open MATLAB and run `./relaxOpt/relaxedOpt.m`.

## License

MIT

