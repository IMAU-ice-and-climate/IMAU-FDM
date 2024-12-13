# Setting up ECMWF to run the IMAU Firn Densification Model

## Resources
[Training Manual](https://confluence.ecmwf.int/display/UDOC/Atos+HPCF+and+ECS+Introduction+Tutorial)
[Monitoring compute time/resources](https://hpc-usage.ecmwf.int/sbu-accounting/my-accounts)

## Pre-requisites
1. request access to IMAU ECMWF account
2. set up multi-factor authentication (login -> account -> change password/OTP)

## Logging in for the first time
[ECMWF how to login documentation](https://www.ecmwf.int/en/computing/access-computing-facilities/how-log)

### SSH

1. install teleport <=13 from [website](https://goteleport.com/download/) (do not use brew! installs version > 13)
2. `teleport-login` (might need to add `--force-clean`)
3. `mkdir ~/.ssh` and `touch ~/.ssh/config`if don't already exist
4. `open ~/.ssh/config` and edit ssh config file, replacing ecmwfusername with your username and user.address@somewhere.com with your email

```
Match host jump.ecmwf.int exec "tsh status --proxy %h >/dev/null 2>&1 || tsh --proxy %h login"
Host jump.ecmwf.int a?-* a??-* hpc-* hpc2020-* ecs-*
  User ecmwfusername
  IdentityFile ~/.tsh/keys/jump.ecmwf.int/user.address@somewhere.com
  CertificateFile ~/.tsh/keys/jump.ecmwf.int/user.address@somewhere.com-ssh/jump.ecmwf.int-cert.pub
  HostKeyAlgorithms +ssh-rsa*,rsa-sha2-512
  PubkeyAcceptedKeyTypes +ssh-rsa*
  ServerAliveInterval 60
  TCPKeepAlive yes
 
Host a?-* a??-* hpc-* hpc2020-* ecs-*
  ProxyJump jump.ecmwf.int
```
*note: may need to append `Host ecs-login hpc-login` below last line*
5. run `ssh hpc-login`
6. once you're logged in, run `ssh-key-setup`

### Now that you're in, set up your account for success

1. Copy and source bash profile from user `rukv` _note: an example bashrc profile will be uploaded to github soon.._
```
cd $HOME
cd ../rumb/
cp .bashrc $HOME/.bashrc
source .bashrc
```
2. Transfer or clone FDM package.
- Transfer can be done using rsync, see section below.
- Cloning can bed one from github, e.g. `git clone https://github.com/brils001/IMAU-FDM`. May need to set up Github SSH or personal token first.
3. Compile `Program/Rundir/readpointlist.f90`
- 

## File editing, transfer, and version control

### Using rsync to transfer files

`rsync -avz mydataset hpc-login:PATH`
e.g. `rsync -avz testfile.txt hpc-login:/ec/res4/scratch/USERNAME`; replace USERNAME with your own

*note:* the file paths in the [tutorial](https://confluence.ecmwf.int/display/UDOC/HPC2020%3A+File+transfers) are wrong, use the paths that show up when you login into the node:
```
[ECMWF-INFO -ecprofile] $SCRATCH=/ec/res4/scratch/USERNAME
[ECMWF-INFO -ecprofile] $PERM=/perm/USERNAME
[ECMWF-INFO -ecprofile] $HPCPERM=/ec/res4/hpcperm/USERNAME
```

### Github SSH/personal token for cloning

#### Personal token
1. To create a new personal access token, go to www.githhub.com, then _account_ → _developer settings_ → _new personal token_. Classic token is fine.
2. Copy and save the token somewhere.
2. Any repo can now be cloned using the token: `git clone https://TOKEN@github.com/USERNAME/REPOSITORY.git`

#### SSH
*coming soon*

### Using rmate to edit files in dekstop sublime
_Instructions to open a tunnel between the HPC and your desktop version of Sublime (a text editor). Some people find this much easier to use than vim._
1. Install [sublime](https://www.sublimetext.com/)
2. In Sublime, install rsub using [the package manager](https://packagecontrol.io/installation): cmd/ctrl+shift+p -> type `Package Control: Install Package` -> type `rsub` and hit enter
3. Restart Sublime
4. In your terminal, open `.ssh/config` and add `RemoteForward 52698 127.0.0.1:52698` under`Host jump.ecmwf.int a?-* a??-* hpc-* hpc2020-* ecs-*`.
5. SSH into ECMWF 
6. Run the following to download `rmate`:
```
curl -Lo ~/bin/rmate https://raw.githubusercontent.com/textmate/rmate/master/bin/rmate
chmod a+x ~/bin/rmate
```
7. Add `chmod a+x ~/bin/rmate` and `export PATH="$PATH:$HOME/bin"` to your `bashrc`
8. Now you can type `rmate FILENAME` and edit the file in Sublime on your system. 

