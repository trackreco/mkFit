# mictest: a repository for vectorized, parallelized charged particle track reconstruction

Intro: Below is a short README on setup steps, code change procedures, and some helpful pointers. Please read this thoroughly before checking out the code! As this is a markdown file, it is best viewed via a web browser.

### Outline
1) Test platforms
2) How to checkout the code
3) How to run the code
4) How to make changes to the main development branch
5) How to run the benchmarking and validation suite
6) Submit an issue
7) Condensed description of code
8) Other helpful README's in the repository
9) CMSSW integration
   1) Considerations for `mkFit` code
   2) Instructions to use `mkFit` from CMSSW
      1) Build `mkFit`
      2) Set up `mkFit` as an external
      3) Pull CMSSW code and build
      4) Use `mkFit` in InitialStep of CMS offline tracking
      5) Producing MultiTrackValidator plots
10) Other useful information
    1) Important Links
    2) Tips and Tricks
       1) Missing Libraries and Debugging
       2) SSH passwordless login for benchmarking and web scripts
    3) Acronyms/Abbreviations

## Section 1: Test platforms

- **phi1.t2.ucsd.edu**: _Intel Xeon Sandy Bridge_ (referred to as SNB, phiphi, phi1), model number: E5-2620 [https://ark.intel.com/products/64594/Intel-Xeon-Processor-E5-2620-15M-Cache-2_00-GHz-7_20-GTs-Intel-QPI]
- **phi2.t2.ucsd.edu**: _Intel Xeon Phi Knights Landing_ (referred to as KNL, phi2), model number: 7210 [https://ark.intel.com/products/94033/Intel-Xeon-Phi-Processor-7210-16GB-1_30-GHz-64-core] 
- **phi3.t2.ucsd.edu**: _Intel Xeon Skylake-SP Gold_ (referred to as SKL-Au, SKL-SP, phi3), model number: 6130 [https://ark.intel.com/products/120492/Intel-Xeon-Gold-6130-Processor-22M-Cache-2_10-GHz]
- **lnx4108.classe.cornell.edu**: _Intel Xeon Skylake-SP Silver_ (referred to as SKL-Ag, SKL-SP, lnx4108), model number: 4116 [https://ark.intel.com/products/120481/Intel-Xeon-Silver-4116-Processor-16_5M-Cache-2_10-GHz] 
- **GPUs**: to be filled out

phi1, phi2, and phi3 are all managed across a virtual login server and therefore the home user spaces are shared. phi1, phi2, phi3, and lnx4108 also have /cvmfs mounted so you can source the environment needed to run the code. 

The main development platform is phi3. This is the recommended machine for beginning development and testing. Login into any of the machines is achieved through ```ssh -X -Y <phi username>@phi<N>.t2.ucsd.edu```. It is recommended that you setup ssh key forwarding on your local machine so as to avoid typing in your password with every login, and more importantly, to avoid typing your password during the benchmarking (see Section 10.ii.b).

## Section 2: How to checkout the code

The master development branch is ```devel``` on https://github.com/cerati/mictest (referred to as ```cerati/devel``` for the remainder of the README). This is a private repository, as are all forks of this repository. Development for mictest is done on separate branches within a forked repository. Since Giuseppe is politely hosting the main repo on his account, make sure to fork the repository to your own account first (using the "Fork" option at the top of the webpage), and push any development branches to your own forked repo first.

Once forked, checkout a local copy by simply doing a git clone:

```
git clone git@github.com:<user>/mictest
```

where ```<user>``` is your GH username if renamed your remote to your username. Otherwise ```<user>``` will be ```origin```.

If you wish to add another user's repo to your local clone, do:

```
git remote add <user> git@github.com:<user>/mictest
```

This is useful if you want to submit changes to another user's branches. To checkout a remote branch, do:

```
git fetch <user>
git fetch <user> <branch>
git checkout -b <branch> <user>/<branch>
```

## Section 3: How to run the code

As already mentioned, the recommended test platform to run the code is phi3. Checkout a local repo on phi3 from your forked repo. To run the code out-of-the-box from the main ```devel``` branch, you will first need to source the environment:

```
source xeon_scripts/init-env.sh
```

You are free to put the lines from this script in your login scripts (.bashrc, .bash_profile, etc). However, encapsulate them within a function and then call that function upon logging into phi3. We want clean shells before launching any tests. Therefore, if you have any setup that sources something, disable it and do a fresh login before running any tests! 

Now compile the code:

```
make -j 32 AVX_512:=1
```

To run the code with some generic options, do:

```
./mkFit/mkFit --cmssw-n2seeds --input-file /data2/slava77/samples/2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/memoryFile.fv3.clean.writeAll.recT.072617.bin --build-ce --num-thr 64 --num-events 20
```

Consult Sections 7-8 for where to find more information on descriptions of the code, which list resources on where to find the full set of options for running the code.

There are ways to run this code locally on macOS. Instructions for how to to this will be provided later. You will need to have XCode installed (through the AppStore), XCode command line tools, a ROOT6 binary (downloaded from the ROOT webpage), as well as TBB (through homebrew). 

## Section 4: How to make changes to the main development branch

Below are some rules and procedures on how to submit changes to the main development branch. Although not strictly enforced through settings on the main repo, please follow the rules below. This ensures we have a full history of the project, as we can trace any changes to compute or physics performance that are introduced (whether intentional or unintentional). 

**Special note**: Do not commit directly to ```cerati/devel```! This has caused issues in the past that made it difficult to track down changes in compute and physics performance. Please always submit a Pull Request first, ensuring it is reviewed and given the green light before hitting "Merge pull request". 

1. Checkout a new branch on your local repo: ```git checkout -b <branch>```
2. Make some changes on your local repo, and commit them to your branch: ```git commit -m "some meaningful text describing the changes"```
3. If you have made multiple commits, see if you can squash them together to make the git history legibile for review. If you do not know what you are doing with this, make sure to save a copy of the local branch as backup by simplying checking out a new branch from the branch you are with something like: ```git checkout -b <branch_copy>```. Git provides a tutorial on squashing commits: [https://git-scm.com/book/en/v2/Git-Tools-Rewriting-History]
4. Ensure you have pulled down the latest changes from the main development branch merged into your local development branch. ```git merge cerati devel``` can make a mess, so the preferred option is ```git rebase --onto <new_base_hash> <old_base_hash> <branch>```. CMSSW provides a nice explanation of this rebase option: [https://cms-sw.github.io/tutorial-resolve-conflicts.html]
5. Test locally!
   1. If you have not done so, clone your forked repo onto phi3, checking out your new branch.
   2. Source the environment for phi3 as explained in Section 3.
   3. Compile test: ```make -j 32 AVX_512:=1```. Fix compilation errors if they are your fault or email the group / person responsible to fix their errors! 
   4. Run benchmark test: ```./mkFit/mkFit --cmssw-n2seeds --input-file /data2/slava77/samples/2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/memoryFile.fv3.clean.writeAll.recT.072617.bin --build-ce --num-thr 64 --num-events 20```. Ensure the test did not crash, and fix any segfaults / run-time errors! 
   5. Compile with ROOT test: ```make -j 32 AVX_512:=1 WITH_ROOT:=1```. Before compiling, make sure to do a ```make distclean```, as we do not want conflicting object definitions. Fix errors if compilation fails.
   6. Run validation test:  ```./mkFit/mkFit --cmssw-n2seeds --input-file /data2/slava77/samples/2017/pass-4874f28/initialStep/PU70HS/10224.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017PU_GenSimFullINPUT+DigiFullPU_2017PU+RecoFullPU_2017PU+HARVESTFullPU_2017PU/a/memoryFile.fv3.clean.writeAll.recT.072617.bin --build-ce --num-thr 64 --num-events 20 --backward-fit-pca --cmssw-val-fhit-bprm```. Ensure the test did not crash! 
6. Run the full benchmarking + validation suite on all platforms: follow procedure in Section 5 (below)! If you notice changes to compute or physics performance, make sure to understand why! Even if you are proposing a technical two-line change, please follow this step as it ensures we have a full history of changes.
7. Prepare a Pull Request (PR)
   1. Push your branch to your forked repo on GitHub: ```git push <forked_repo_name> <branch>```
   2. Navigate to the main GH: https://github.com/cerati/mictest
   3. Click on "New Pull Request"
   4. Click on "Compare across forks", and navigate to your fork + branch you wish to merge as the "head fork + compare"
   5. Provide a decent title, give a brief description of the proposed commits. Include a link to the benchmarking and validation plots in the description. If there are changes to the compute or physics performance, provide an explanation for why! If no changes are expected and none are seen, make sure to mention it.
   6. (Optional) Nominate reviewers to check over the proposed changes.
   7. Follow up on review comments! After pushing new commits to your branch, repeat big steps 5 and 6 (i.e. test locally and re-run the validation). Post a comment to the PR with the new plots.
   8. Once given the green light, you can hit "Merge Pull Request", or ask someone else to do it.

## Section 5: How to run the benchmark and validation suite

Notes on nomenclature:
- "benchmark": these are the compute performance tests (i.e. time and speedup)
- "validation": these are the physics performance tests (i.e. track-finding efficiency, fake rate, etc.)

We often use these words interchangibly to refer to the set of benchmark and validation tests as a single suite. So if you are asked to "run the benchmarking" or "run the validation": please run the full suite (unless specifically stated to run one or the other). In fact, the main scripts that run the full suite use "benchmark" in their name, even though they may refer to both the running of the compute and physics performance tests and plot comparisons.

When submitting a PR or preparing for a conference, please run the full suite from phi3 with a clean login: make sure nothing has been sourced to set up the environment. The scripts will set the environment as needed through ```./xeon_scripts/common-variables.sh ${suite}```. 

The main script for running the full suite can be launched with:

```
./xeon_scripts/runBenchmark.sh ${suite}
```

There are three options for running the full suite by passing one of the three strings to the parameter ```${suite}```:
- ```full``` : runs compute and physics tests for all track finding routines (BH, STD, CE, FV)
- ```forPR``` : runs compute and physics tests for track finding routines used for comparisons in pull requests (default setting: BH and CE for benchmarks, STD and CE for validation)
- ```forConf``` : runs compute and physics tests for track finding routines used for conferences only (currently only CE)

The ```full``` option currently takes little over an hour to run.  Make sure the machines are quiet before launching any tests: we don't want to disturb someone who already is testing! 

Inside the main script, tests are submitted for phi1, phi2, and phi3 concurrently, by tarring up the local repo, sending the tarball to a disk space on the remote platform, compiling the untarred directory natively on the remote platform, and then sending back the results to be collected on phi3. It should be noted that the tests for phi3 are simply run on in the user home directory when logged into phi3 (although we could in principle ship the code to the work space disk on phi3). Because we run the tests for phi3 in the home directory, which is shared by all three machines, we pack and send the code to a remote _disk_ space _before_ launching the tests on phi3 from the home directory. The scripts that handle the remote testing are: 

```
./xeon_scripts/tarAndSendToRemote.sh ${remote_arch} ${suite}
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build-remote.sh ${ben_arch} ${suite}
```

When these scripts are called separately to run a test on particular platform, one of three options must be specified for ```${remote_arch}``` or ```${ben_arch}```: ```SNB```, ```KNL```, or ```SKL-SP```. The main script ```runBenchmark.sh``` will do this automatically for all three platforms. If the code has already been compiled on a given machine, it is sufficient to run:

```
./xeon_scripts/benchmark-cmssw-ttbar-fulldet-build.sh ${ben_arch} ${suite}
```

The appropriate strings should appear in place of ```${ben_arch}``` and ```${suite}```.

Within the main ```./xeon_scripts/runBenchmark.sh``` script, there are two other scripts that make performance plots from the log files of compute performance tests:

```
./plotting/benchmarkPlots.sh ${suite} 
./plotting/textDumpPlots.sh ${suite}
```

The first will produce the time and speedup plots, while the second produces distributions of basic kinematic quantites of the candidate track collections, comparing the results across the different platforms and different number of vector units and threads. Ideally, in plots from the second script, the distributions should have all points lie on top of each other: there should be no dependency on platform or parallelization/vectorization setting for a specific track-finding routine.

The main physics performance script that is run is:

```
./val_scripts/validation-cmssw-benchmarks.sh ${suite}
```

This script will run the validation on the building tests specified by the ```${suite}``` option. It will also produce the full set of physics performance plots and text files detailing the various physics rates.

It should be mentioned that each of these scripts within ```./xeon_scripts/runBenchmark.sh``` can be launched on their own, as again, they each set the environment and run tests and/or plot making. However, for simplicity's sake, it is easiest when prepping for a PR to just run the master ```./xeon_scripts/runBenchmark.sh```.  If you want to test locally, it is of course possible to launch the scripts one at a time.

Note that you need to have SSH-forwarding set up to avoid having to type your password every time ```./xeon_scripts/runBenchmark.sh``` needs to copy files back and forth.

After running the full suite, there is an additional set of scripts within the ```web/``` directory for organizing the output plots and text files for viewing them on the web.  The main script is:

```
./web/move-benchmarks.sh ${outdir_name} ${suite} ${afs_or_eos}
```

where again, ```${suite}``` defaults to ```forPR```. ```${outdir_name}``` will be the top-level directory where the output is collected and eventually shipped to LXPLUS. This will call ```./web/collectBenchmarks.sh ${outdir_name} ${suite}```, which will sort the files, and then ```./web/tarAndSendToLXPLUS.sh ${outdir_name} ${suite} ${afs_or_eos}```, which  packs up the top-level output dir and copies it to either an /afs or /eos userspace on LXPLUS. This will also run another script remotely to copy ```web/index.php``` into each directory to have a nice web GUI for the plots. Make sure to read the ```web/README_WEBPLOTS.md``` first to setup an /afs or /eos web directory. 

The option ```${afs_or_eos}``` takes either of the following arguments: ```afs``` or ```eos```, and defaults to ```afs```. The mapping of usernames to /afs or /eos spaces is in ```./xeon_scripts/common-variables.sh```. If an incorrect string is passed, the script will exit. 

**IMPORTANT DISCLAIMER**

There is a script: ```./xeon_scripts/trashSKL-SP.sh``` that is run at the very end of the ```./web/move-benchmarks.sh``` script that will delete log files, pngs, validation directories, root files, and the neat directory created to house all the validation plots.  This means that if the scp fails, the plots will still be deleted locally, and you will be forced to re-run the whole suite!!  You can of course comment this script out if this bothers you.

**FINAL WORD ON AFS**

AFS is being phased out at CERN. The web scripts currently accomodate both options. However, if you have your own website you want to post plots to, then feel free to use ```./web/collect-benchmarks.sh``` to tidy up the plots into neat directories, and then send them where you want.

## Section 6: Submit an issue

It may so happen that you discover a bug or that there is a known problem that needs further discussion outside of private emails/the main list-serv. If so, make sure to open issue on the main repo by clicking on "Issues" on GH, then "Open an issue".  Provide a descriptive title and a description of the issue. Provide reference numbers to relevant PRs and other Issues with"#<number>".  Include a minimal working example to reproduce the problem, attaching log files of error messages and/or plots demonstrating the problem. 

Assign who you think is responsible for the code (which could be yourself!). If you have an idea that could solve the problem: propose it! If it requires a large change to the code, or may hamper performance in either physics or computing, make sure to detail the pros and cons of different approaches. 

Close an issue after it has been resolved, providing a meaningful message + refence to where/how it was resolved.

## Section 7: Condensed description of code

### mkFit/mkFit.cc

This file is where the ```main()``` function is called for running the executable ```./mkFit/mkFit```. The ```main()``` call simply setups the command line options (and lists them), while the meat of the code is called via ```test_standard()```. Some of the command line options will set global variables within mkFit.cc, while others will set the value of variables in the ```Config``` nampespace. Options that require strings are mapped to via enums in the code, with the mapping specified via global functions at the top of mkFit.cc

```test_standard()``` does the majority of the work: running the toy simulation, reading or writing binary files, and running the various tests. The outer loop is a TBB parallel-for over the number of threads used for running multiple-events-in-flight (MEIF). The default is one event in flight. The inner loop is over the number of events specified for that thread. The number of events in total to run over can be specified as a command line option. When running multiple-events-in-flight, in order to have reasonable statistics from variable load from different events, it is advised to have at least 20 events per thread.  When we refer to "total loop time" of the code, we are timing the inner loop section for each event, which includes I/O. However, for the sake of the plots, we simply sum the time for all events and all threads, and divide by the number of events run to obtain an average per event time.

Within the inner loop, a file is read in, then the various building and fitting tests are run. At the end of each event there is optional printout, as well as at the end of all tthe events for a thread. If running the validation with multiple-events-in-flight is enabled, you will have to ```hadd``` these files into one file before making plots. This is handled automatically within the scripts. 

### mkFit/buildtestMPlex.[h,cc]

This code calls the various building routines, setting up the event, etc. The functions defined here are called in mkFit.cc. Functions called within this file are from MkBuilder.

### mkFit/MkBase.h + mkFit/MkFitter.[h,cc] + mkFit/MkFinder.[h,cc]

MkFinder and MkFitter derive from MkBase. High-level code for objects used by building and fitting routines in mkFit. These objects specify I/O operations from standard format to Matriplex format for different templated Matriplex objects (see Matrix[.h,.cc] for template definitions). 

### mkFit/MkBuilder.[h,cc]

Specifies building routines, seed prepping, validation prepping, etc. Code for building and backward fit routines using MkFinders, while seed fitting uses MkFitters. Objects from Event object are converted to their Matriplex-ready formats. Uses the layer plan to navigate which layer to go to for each track. Foos for the navigation are defined in SteerinParams.h.

### Math/ directory

Contains SMatrix headers, used for some operations on track objects (mostly validation and deprecated SMatrix building code -- see below).

### Matriplex/ directory

Contains low-level Matriplex library code for reading/writing into matriplex objects as well as elementary math operations (add, multiply). Includes perl scripts to autogenerate code based on matrix dimension size.

### Geoms/ dir + TrackerInfo.[h,cc]

Geometry plugin info. TrackerInfo setups classes for layer objects. Geoms/ dir contains the actual layout (number scheme, layer attributes, etc) for each of the different geoemetries.

### mkFit/PropagationMPlex.[h,cc,icc] + mkFit/KalmanUtilsMPlex.[h,cc,icc]

Underlying code for propagation and Kalman upate (gain) calculations in Matriplex form. The .icc files contain the low-level computations. Chi2 computations specified in KalmanUtilsMPlex.

### mkFit/CandCloner.[h,cc]

Code used in Clone Engine for bookkeeping + copying candidates after each layer during building. 

### mkFit.HitStructures.[h,cc]

Specifies MkBuilder + Matriplex friendly data formats for hits. Hits are placed in these containers before building.

### Event.[h,cc]

Most of the code is vestigial (see below). However, the Event object is a container for the different track collections and hit collection. There is code for seed processing, namely cleaning. There is also code relevant for validation and validation prep for different track collections.

### Hit.[h,cc] + Track.[h,cc]

Contain the Hit, Track, and TrackExtra classes. These are the "native" formats read from the binary file (read in from the Tracking NTuple). In principle, since we are planning to migrate to CMSSW eventually, these classes (as well Event) may be trimmed to just read straight from CMSSW native formats.

- Hit object contains hit parameters, covariance, and a global ID. The global ID is used for gaining more information on the MC generation of that hit.
- Track object is simply the track parameters, covariance, charge, track ID, and hit indices + layers. 
- TrackExtra contains additional information about each track, e.g. associated MC info, seed hits, etc. A Track's TrackExtra is accessed through the track label, which is the index inside the vector of tracks. 

### Config.[h,cc]

Contains the Config namespace. Specifies configurable parameters for the code. For example: number of candidates to create for each track, chi2 cut, number of seeds to process per thread, etc. Also contains functions used for dynamically setting other parameters based on options selected. 

Tracker Geometry plugin also initialized here.

### Validation code

Described in validation manifesto. See Section 8 for more info on manifesto.

### TO DO

- flesh out sections as needed
- GPU specific code

### Vestigial code

There are some sections of code that are not in use anymore and/or are not regularly updated. A short list is here:
- main.cc : Old SMatrix implementation of the code, which is sometimes referred to as the "serial" version of the code.
- USolids/ : Directory for implementing USolids geometry package. Originally implemented in SMatrix code.
- seedtest[.h,.cc] : SMatrix seeding
- buildtest[.h,.cc] : SMatrix building
- fittest[.h,.cc] : SMatrix fitting
- ConformalUtils[.h,.cc] : SMatrix conformal fitter for seeding/fitting
- (possibly) Propagation[.h,.cc] : currently in use by the currently defunct Simulation[.h,.cc]. In reality, will probably move simulation code to MPlex format, which will deprecate this code.
- KalmanUtils[.h,.cc] : SMatrix Kalman Update
- mkFit/seedtestMPlex[.h,.cc] and all code in MkBuilder[.h,.cc] related to finding seeds with our own algorithm
- mkFit/ConformalUtils[.h,.cc] : used by the seeding, although could be revived for fitting
- additional val_scripts/ and web/ scripts not automatically updated outside of main benchmarking code
- mtorture test/ code 

## Section 8: Other helpful README's in the repository

Given that this is a living repository, the comments in the code may not always be enough. Here are some useful other README's within this repo:
- afer compiling the code, do: ```./mkFit/mkFit --help``` : Describes the full list of command line options, inputs, and defaults when running mkFit. The list can also be seen in the code in mkFit/mkFit.cc, although the defaults are hidden behind Config.[h,cc], as well as mkFit.cc.
- cmssw-trackerinfo-desc.txt : Describes the structure of the CMS Phase-I geometry as represented within this repo.
- index-desc.txt : Desribes the various hit and track indices used by different sets of tracks throughout the different stages of the read in, seeding, building, fitting, and validation.
- validation-desc.txt : The validation manifesto: (somewhat) up-to-date description of the full physics validation suite. It is complemented by a somewhat out-of-date code flow diagram, found here: https://indico.cern.ch/event/656884/contributions/2676532/attachments/1513662/2363067/validation_flow_diagram-v4.pdf
- web/README_WEBPLOTS.md : A short markdown file on how to setup a website with an AFS or EOS directory on LXPLUS (best when viewed from a web browser, like this README).

## Section 9: CMSSW integration

The supported CMSSW version is currently `10_2_0_pre3`. The
integration of `mkFit` in CMSSW is based on setting it up as a CMSSW
external.

### Section 9.i: Considerations for `mkFit` code

The multi-threaded CMSSW framework, and the iterative nature of CMS
tracking impose some constraints on `mkFit` code (that are not all met
yet). Note that not all are mandatory per se, but they would make the
life easier for everybody.

* A single instance of `mkFit` should correspond to a single track building iteration
* There should be no global non-const variables
  - Currently there are non-const global variables e.g. in `Config` namespace
* All iteration-specific parameters should be passed from CMSSW to `mkFit` at run time

### Section 9.ii: Instructions to use `mkFit` from CMSSW

#### Section 9.ii.a: Build `mkFit`

To be used from CMSSW the `mkFit` must be built with the CMSSW
toolchain. Assuming you are in an empty directory, the following
recipe will set up a CMSSW developer area and a `mictest` area there,
and compile `mkFit` using the CMSSW toolchain.

```bash
cmsrel CMSSW_10_2_0_pre3
pushd CMSSW_10_2_0_pre3/src
cmsenv
git cms-init
scram setup icc-ccompiler
source $(scram tool tag icc-ccompiler ICC_CCOMPILER_BASE)/bin/iccvars.sh intel64
popd
git clone git@github.com:cerati/mictest
pushd mictest
TBB_PREFIX=$(dirname $(cd $CMSSW_BASE && scram tool tag tbb INCLUDE)) make -j 12
popd
```

#### Section 9.ii.b: Set up `mkFit` as an external

Assuming you are in the aforementioned parent directory, the following
recipe will create a scram tool file, and set up scram to use it

```bash
pushd CMSSW_10_2_0_pre3/src
cat <<EOF >mkfit.xml
<tool name="mkfit" version="1.0">
  <client>
    <environment name="MKFITBASE" default="$PWD/../../mictest"/>
    <environment name="LIBDIR" default="\$MKFITBASE/lib"/>
    <environment name="INCLUDE" default="\$MKFITBASE"/>
  </client>
  <runtime name="MKFIT_BASE" value="\$MKFITBASE"/>
  <lib name="MicCore"/>
  <lib name="MkFit"/>
</tool>
EOF
scram setup mkfit.xml
cmsenv
```

#### Section 9.ii.c: Pull CMSSW code and build

The following recipe will pull the necessary CMSSW-side code and build it

```bash
# in CMSSW_10_2_0_pre3/src
git cms-remote add makortel
git fetch makortel
git checkout -b mkfit_1020pre3 makortel/mkfit_1020pre3
git cms-addpkg $(git diff $CMSSW_VERSION --name-only | cut -d/ -f-2 | uniq)
git cms-checkdeps -a
scram b -j 12
```

#### Section 9.ii.d: Use `mkFit` in InitialStep of CMS offline tracking

The example below uses 2018 tracking-only workflow

```bash
# Generate configuration
runTheMatrix.py -l 10824.1 --apply 2 --command "--customise RecoTracker/MkFit/customizeInitialStepToMkFit.customizeInitialStepToMkFit" -j 0
cd 10824.1*
# edit step3*RECO*.py to contain your desired (2018 RelVal MC) input files
cmsRun step3*RECO*.py
```

The customize function replaces the initialStep track building module
with `mkFit`. In principle the customize function should work with any
reconstruction configuration file.

By default `mkFit` is configured to use Clone Engine with N^2 seed
cleaning, and to do the backward fit (to the innermost hit) within `mkFit`.

#### Section 9.ii.e: Producing MultiTrackValidator plots

The `step3` above runs also the [MultiTrackValidator](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMultiTrackValidator).

To produce the plots, first run the DQM harvesting step

```bash
cmsRun step4_HARVESTING.py
```

which produces a `DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root` file that contains all the histograms. Rename the file to something reflecting the contents, and run

```bash
makeTrackValidationPlots.py --extended --limit-tracking-algo initialStep <DQM file> [<another DQM file> ...]
```

The script produces a directory `plots` that can be copied to any web area.


## Section 10: Other useful information

### Section 10.i: Important Links
- Main development GitHub: https://github.com/cerati/mictest
- Our project website: https://trackreco.github.io
- Out-of-date and longer used twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MicTrkRnD
- Indico meeting page: https://indico.cern.ch/category/8433
- Vidyo room: Parallel_Kalman_Filter_Tracking
- Email list-serv: mic-trk-rd@cern.ch

### Section 10.ii: Tips and Tricks

#### Section 10.ii.a: Missing Libraries and Debugging

When sourcing the environment on phi3 via ```source xeon_scripts/init-env.sh```, some paths will be unset and access to local binaries may be lost. For example, since we source ROOT (and its many dependencies) over CVMFS, there may be some conflicts in loading some applications. In fact, the shell may complain about missing environment variables (emacs loves to complain about TIFF). The best way around this is to simply use CVMFS as a crutch to load in what you need.

This is particularly noticeable when trying to run a debugger. To compile the code, at a minimum, we must source icc + toolkits that give us libraries for c++14. We achieve this through the dependency loading of ROOT through CVMFS (previously, we sourced devtoolset-N to grab c++14 libraries). 

After sourcing and compiling and then running only to find out there is some crash, when trying to load ```mkFit``` into ``gdb`` via ```gdb ./mkFit/mkFit```, it gives rather opaque error messages about missing Python paths.

This can be overcome by loading ```gdb``` over CVMFS: ```source /cvmfs/cms.cern.ch/slc7_amd64_gcc630/external/gdb/7.12.1-omkpbe2/etc/profile.d/init.sh```. At this point, the application will run normally and debugging can commence.

#### Section 10.ii.b: SSH passwordless login for benchmarking scripts and web scripts

When running the benchmarks, a tarball of the working directory will be ```scp```'ed to phi2 and phi1 before running tests on phi3. After the tests complete on each platform, the log files will be ```scp```'ed back to phi3 concurrently. If you do not forward your ssh keys upon login to phi3, you will have to enter your password when first shipping the code over to phi2 and phi1, and also, at some undetermined point, enter it again to receive the logs.

With your favorite text editor, enter the text below into ```~/.ssh/config``` on your local machine to avoid having to type in your password for login to any phi machine (N.B. some lines are optional):

```
Host phi*.t2.ucsd.edu
     User <phi* username>
     ForwardAgent yes
# lines below are for using X11 on phi* to look at plots, open new windows for emacs, etc.
     ForwardX11 yes
     XAuthLocation /opt/X11/bin/xauth
# lines below are specific to macOS	     
     AddKeysToAgent yes 
     UseKeychain yes
```

After the benchmarks run, you may elect to use the ```web/``` scripts to transfer plots to CERN website hosted on either LXPLUS EOS or AFS. The plots will be put into a tarball, ```scp```'ed over, and then untarred remotely via ```ssh```. To avoid typing in your password for the ```web/``` scripts, you will need to use a Kerberos ticket and also modify your ```.ssh/config``` file in your home directory on the _phi_ machines with the text below:

```
Host lxplus*.cern.ch
     User <lxplus username>
     ForwardAgent yes
     ForwardX11 yes
     GSSAPIAuthentication yes
     GSSAPIDelegateCredentials yes
```

The last two lines are specific to Kerberos's handling of ssh, which is installed on all of the _phi_ machines. In order to open a Kerberos ticket, you will need to do:

```
kinit -f <lxplus username>@CERN.CH
```

and then enter your LXPLUS password. Kerberos will keep your ticket open for a few days to allow passwordless ```ssh``` into LXPLUS. After the ticket expires, you will need to enter that same command again. So, even if you only send plots once every month to LXPLUS, this reduces the number of times of typing in your LXPLUS password from two to one :).

### Section 10.iii: Acronyms/Abbreviations:
- AVX: Advanced Vector Extensions [flavors of AVX: AVX, AVX2, AVX512]
- BH: Best Hit
- CCC: Charge Cluster Cut
- CE: Clone Engine
- CMS: Compact Muon Solenoid
- CMSSW: CMS Software
- CMSSWVal: CMSSWTrack Validation, use cmssw tracks as reference set of tracks for association
- FV: Full Vector
- GH: GitHub
- GUI: Graphical User Interface
- MEIF: Multiple-Events-In-Flight
- MP: Multi-Processing
- N^2: Local seed cleaning algorithm developed by Mario and Slava
- PR: Pull Request
- Reco: Reconstruction
- SimVal: SimTrack validation, use simtracks as reference set of tracks for association
- STD: Standard
- TBB: (Intel) Threaded Building Blocks
- TH: Threads
- VU: (loosely) Vector Units

Glossary of acronyms from CMS: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGlossary
