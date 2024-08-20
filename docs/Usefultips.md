# Useful tips and commands

Here are some useful command lines you may need to work on ENS local computer and the IFB-core Cluster during the practical.

You can find more detailed information on [IFB-core Cluster Documentation](https://ifb-elixirfr.gitlab.io/cluster/doc/).

***
## How to download/upload your data from the cluster

To download files from the cluster to your current directory
```
scp  '<your login>@core.cluster.france-bioinformatique.fr:/<absolute path to your file>' .
```
   
To download a folder from the cluster to your current directory
```
scp -r  '<your login>@core.cluster.france-bioinformatique.fr:/<absolute path to your folder>' .
```
    
To upload a file to the cluster
```
scp '<path to your local file>' '<your login>@core.cluster.france-bioinformatique.fr:/<absolute path to the target folder>'
```

***
## IFB-core cluster connexion

Beware of passwords that contains a '^' sign, you will need to copy/paste the password from a text file in order to work on ENS computers.

If you get an error message because your login is too long  'too long for Unix domain socket' you could replace the server name by its IP address: 192.54.201.181
```
ssh <your login>@192.54.201.181
```
