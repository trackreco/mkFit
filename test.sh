source xeon_scripts/common-variables.sh ${suite} ${useLNX} ${lnxuser}
source xeon_scripts/init-env.sh

HOST=${LNXG_HOST}

check_settings=false
echo "--------Showing System Settings--------"
# unzip tarball remotely
echo "Untarring repo on ${remote_arch} remotely"
SSHO ${HOST} bash -c "'
echo "--------Showing System Settings--------"
##### Check Settings #####
echo "turbo status: "$(cat /sys/devices/system/cpu/intel_pstate/no_turbo)
echo "scaling governor setting: "$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)
echo ${check_settings}
echo "--------End System Settings ------------"

echo ${check_settings}
if [ ${check_settings} ] ;
then
echo "Ensuring correct settings"
if [[ $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor) != "performance" ]]
then
echo "performance mode is OFF. Exiting"
exit 1
fi
if [[ $(cat /sys/devices/system/cpu/intel_pstate/no_turbo) == "0" ]]
then
echo "Turbo is ON. Exiting"
exit 1
fi
fi
sleep 3 ## so you can see the settings
'"
bad=$(SSHO ${HOST} echo $?)

wait
echo ${bad}
if [ ${bad} -eq 0 ]; then
echo "not killed"
else
echo "killed"
fi
