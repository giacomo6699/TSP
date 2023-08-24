EXE="./tsp"

for SEED in {0..74}; do
    NAME=greedy10_twoopt_1500nodes/tspgreedy10twooptn1500seed${SEED}
	./tsp -seed $SEED -n 1500 -grasp 1 -patch_heur 1 -tl 240 > $NAME.log
done
for SEED in {0..74}; do
    NAME=vns_1500nodes/tspvnsn1500seed${SEED}
	./tsp -seed $SEED -n 1500 -grasp 2 -patch_heur 1 -tl 240 > $NAME.log
done