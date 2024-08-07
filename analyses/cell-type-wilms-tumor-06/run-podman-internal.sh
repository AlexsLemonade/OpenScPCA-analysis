# ids defined in image for the rstudio user, if not define as such, it is not possible to login to RStudio
uid=1000
gid=1000
# subid ranges on host
subuidSize=$(( $(podman info --format "{{ range .Host.IDMappings.UIDMap }}+{{.Size }}{{end }}" ) - 1 ))
subgidSize=$(( $(podman info --format "{{ range .Host.IDMappings.GIDMap }}+{{.Size }}{{end }}" ) - 1 ))

# go to the OpenScPCA-analysis folder before mounting $PWD to keep access to OpenScPCA-analysis/data in the R Session
cd ../.. \

podman run \
  --name maudp_ScPCA_wilms \
  -e RUNROOTLESS=false \
  --uidmap $uid:0:1 --uidmap 0:1:$uid --uidmap $(($uid+1)):$(($uid+1)):$(($subuidSize-$uid)) \
  --gidmap $gid:0:1 --gidmap 0:1:$gid --gidmap $(($gid+1)):$(($gid+1)):$(($subgidSize-$gid)) \
  --group-add=keep-groups \
  --volume=$PWD/:/home/rstudio \
  -e PASSWORD=$PASSWORD \
  -p 8787:8787 \
  -e TZ=Europe/Vienna \
  openscpca/cell-type-wilms-tumor-06:latest
