### Feigenbaum
feigenbaum is a set of project scripts purposed to ensure that our 2 gene ratio experiment actually terminates; It will be composed of projects much like before with a few tweaks.
Specifically,the experiments have different population sizes, which will be done until I find a 'sweet' spot where the simulation equilibriates reliably. Maybe do some tests to come up with a score cut off instead of a generational cut off.

For now, there will be no weight given to expression error.

## STRUCTURE

bin/ contains the higher order scripts that manipulate the control model and execute the output model project scripts.
iterations/ contains the experiment folders for each iteration, and the csv output for each iteration.
paper/ contains what I have completed so for of the paper.
plotting/ contains plot utilities (written in R) and plots.

## INSTRUCTIONS:

		1) Run bin/replace_script.sh and input iteration number (this creates a directory in the iteration/date and iteration/experiments folders with
		the specified iteration name) EX: 1st, 2nd, 3rd, etc
		2) Ensure that the codons are correct in data/model_codons.csv - the script bin/make_experiments.sh changes this file by incrementing 'TAT'
		3) Run bin/make_experiments.sh to create the experiment folders in iterations/experiments
		4) Run sed -e 's/[old name]/[new name]/g' bin/execute_[old name].sh > execute_[new name].sh to make the new execution script
		5) chmod +x bin/execute_[new name].sh to run all experiments.
## EXPERIMENTS:
		1)  experiment has W=0, N=50,  terminated at approx. 500 for 1ratio
		2)  experiment has W=0, N=50,  terminated at approx. 506 for 1ratio
		3)  experiment has W=0, N=50,  terminated at approx. 504 for 1ratio
		4)  experiment has W=0, N=10,  terminated at approx. 35 for 1ratio
		5)  experiment has W=0, N=10,  terminated at approx. 33 for 1ratio
		6)  experiment has W=0, N=50,  terminated at approx. 506 for 1ratio
		7)  experiment has W=0, N=50,  terminated at approx. 505 for 1ratio
		8)  experiment has W=0, N=100, terminated at approx. 505 for 1ratio
		9)  experiment has W=0, N=100, terminated at approx. 504 for 1ratio
		10) experiment has W=0, N=100, terminated at approx. 506 for 1ratio
		11) experiment has W=0, N=100, terminated at approx. 507 for 1ratio
		12) experiment has W=0, N=10,  terminated at approx. 41 for 1ratio
		13) experiment has W=0, N=10,  terminated at approx. 34 for 1ratio
		14) experiment has W=0, N=10,  terminated at approx. 36 for 1ratio

Next, choose parameter N=50 and try multiple weights.

		15) experiment has W=.25, N=50,
		16) experiment has W=.25, N=50,
		17) experiment has W=.25, N=50,
		18) experiment has W=.25, N=50,
		19) experiment has W=.25, N=50,
		20) experiment has W=.50, N=50,
		21) experiment has W=.50, N=50,
		22) experiment has W=.50, N=50,
		23) experiment has W=.50, N=50,
		24) experiment has W=.50, N=50,
		25) experiment has W=.75, N=50,
		26) experiment has W=.75, N=50,
		27) experiment has W=.75, N=50,
		28) experiment has W=.75, N=50,
		29) experiment has W=.75, N=50,
		30) experiment has W=1.0, N=50,
		31) experiment has W=1.0, N=50,
		32) experiment has W=1.0, N=50,
		33) experiment has W=1.0, N=50,
		34) experiment has W=1.0, N=50,

		35) experiment has W=1.25, N=50,
		36) experiment has W=1.25, N=50,
		37) experiment has W=1.25, N=50,
		38) experiment has W=1.25, N=50,
		39) experiment has W=1.25, N=50,
		40) experiment has W=1.5, N=50,
		41) experiment has W=1.5, N=50,
		42) experiment has W=1.5, N=50,
		43) experiment has W=1.5, N=50,
		44) experiment has W=1.5, N=50,
		45) experiment has W=1.75, N=50,
		46) experiment has W=1.75, N=50,
		47) experiment has W=1.75, N=50,
		48) experiment has W=1.75, N=50,
		49) experiment has W=1.75, N=50,
		50) experiment has W=2.0, N=50,
		51) experiment has W=2.0, N=50,
		52) experiment has W=2.0, N=50,
		53) experiment has W=2.0, N=50,
		54) experiment has W=2.0, N=50,
		

Directories have structure:
── experiments 
		└──── pp_amount
			└── pp_rate
					├── 50_expression_error
					│   ├── weight_0
					│   │   ├── 1st_iteration
					│   │   ├── 2nd_iteration
					│   │   ├── 3rd_iteration
					│   │   ├── 4th_iteration
					│   │   └── 5th_iteration
					│   ├── weight_1
					│   │   ├── 1st_iteration
					│   │   ├── 2nd_iteration
					│   │   ├── 3rd_iteration
					│   │   ├── 4th_iteration
					│   │   └── 5th_iteration
					│   ├── weight_2
					│   │   ├── 1st_iteration
					│   │   ├── 2nd_iteration
					│   │   ├── 3rd_iteration
					│   │   ├── 4th_iteration
					│   │   └── 5th_iteration
					│   ├── weight_3
					│   │   ├── 1st_iteration
					│   │   ├── 2nd_iteration
					│   │   ├── 3rd_iteration
					│   │   ├── 4th_iteration
					│   │   └── 5th_iteration
					│   └── weight_99
					│       ├── 1st_iteration
					│       ├── 2nd_iteration
					│       ├── 3rd_iteration
					│       ├── 4th_iteration
					│       └── 5th_iteration
					└── init_expression_error
							└── weight_0
									├── 1st_iteration
									├── 2nd_iteration
									├── 3rd_iteration
									├── 4th_iteration
									└── 5th_iteration


22 directories
