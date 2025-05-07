N_EV=$1
N_PROC=$2

N_EV_PER_PROC=$(($N_EV/$N_PROC))

echo "Number of events per process: ${N_EV_PER_PROC}"
echo "Number of processes: ${N_PROC}"

# generate commands
rm cmds
for i in $(seq 0 ${N_PROC}); do
  j=$(($i+1))
  echo "./toyMCebye ${N_EV_PER_PROC} ${i} ${j}" >> cmds
done

# compile main program
g++ toyMCebye.cxx -o toyMCebye `root-config --cflags --glibs`

# run processes
echo "Running processes..."
cat cmds | xargs -P ${N_PROC} -I CMD bash -c CMD
