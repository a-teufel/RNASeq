#!/bin/bash
#Execute all of the experiments from their relative directories
pwd
(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/100ratio/; bin/execScript.sh)  &
#(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/90ratio/; bin/execScript.sh) &
(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/75ratio/; bin/execScript.sh) &
(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/50ratio/; bin/execScript.sh) &
(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/25ratio/; bin/execScript.sh) &
#(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/10ratio/; bin/execScript.sh) &
#(cd iterations/experiments/pp_rate/50_expression_error/weight_1/3rd_iteration/1ratio; bin/execScript.sh) & wait
