#!/bin/bash


# parse config parameters:
source bash/parse_yaml.sh
eval $(parse_yaml config.yaml CONF_)

# ids defined in image for the rstudio user
uid=1000
gid=1000
# subid ranges on host
subuidSize=$(( $(podman info --format "{{ range .Host.IDMappings.UIDMap }}+{{.Size }}{{end }}" ) - 1 ))
subgidSize=$(( $(podman info --format "{{ range .Host.IDMappings.GIDMap }}+{{.Size }}{{end }}" ) - 1 ))


podman run -d --rm \
  --name ${CONF_project_name}_${USER} \
  -e RUNROOTLESS=false \
  --uidmap $uid:0:1 --uidmap 0:1:$uid --uidmap $(($uid+1)):$(($uid+1)):$(($subuidSize-$uid)) \
  --gidmap $gid:0:1 --gidmap 0:1:$gid --gidmap $(($gid+1)):$(($gid+1)):$(($subgidSize-$gid)) \
  --group-add=keep-groups \
  -p 8080:8787 \
  -e PASSWORD=wordpass \
  -e TZ=Europe/Vienna \
  --volume=$(realpath ${CONF_project_root_host%%*##*( )}):${CONF_project_root} \
  --volume=$(realpath ${CONF_resource_root_host%%*##*( )}):${CONF_resource_root}:ro \
  --volume=$(realpath ${CONF_out_root_host%%*##*( )}):${CONF_out_root} \
  --volume=$(realpath ${CONF_data_root_host%%*##*( )}):${CONF_data_root} \
  ${CONF_project_docker}