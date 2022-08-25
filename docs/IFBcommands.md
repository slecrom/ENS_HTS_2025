# Useful commands to work on the IFB-core Cluster

Here are some useful command lines you may need to work on the IFB-core Cluster during the practical.

You can find more detailed information and video tutorial on [IFB-core Cluster Documentation](https://ifb-elixirfr.gitlab.io/cluster/doc/).

***
## How to download/upload your data from the cluster

To download file/files from the cluster to your current directory
```
scp  '<your login>@core.cluster.france-bioinformatique.fr:/<absolute path to your file>' .
```
   
To download a folder from the cluster to your current directory
```
scp -r  '<your login>@core.cluster.france-bioinformatique.fr:/<absolute path to your folder>' .
```
    
To upload a file to the cluster
```
scp '<path to your local folder>' '<your login>@core.cluster.france-bioinformatique.fr:/<absolute path to the target folder>'
```

***
## How to get information on your current job

```
squeue -u <your login>
```

***
## How to list all your running/pending jobs

```
squeue -u <your login> -t RUNNING
squeue -u <your login> -t PENDING
```

***
## How to cancel jobs

Cancel a specific job
```
scancel <jobid>
```

Cancel all your jobs
```
scancel -u <your login>
```
   
Cancel all your pending jobs
```
scancel -t PENDING -u <your login>
```
