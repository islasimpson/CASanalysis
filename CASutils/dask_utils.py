import dask
import socket

def get_dask_cluster(n_workers, walltime):

    from dask_jobqueue import PBSCluster
    from dask.distributed import Client
    dask.config.set({"distributed.scheduler.worker_saturation":1.0})
    dask.config.set({"optimization.fuse.active": False})
    dask.config.set({
        "distributed.worker.memory.target": 0.6,
        "distributed.worker.memory.spill": 0.7,
        "distributed.worker.memory.pause": 0.8,
        "distributed.worker.memory.terminate": 0.95,
    })
    
    cluster = PBSCluster(
        cores = 1,
        memory = '30GB',
        processes = 1,
        queue = 'casper',
        local_directory = '/glade/derecho/scratch/islas/dask_tmp/',
        resource_spec = 'select=1:ncpus=1:mem=30GB',
        project='P04010022',
        walltime=walltime,
        interface='mgt')
    
    # scale up
    cluster.scale(n_workers)
    #cluster.adapt(minimum=1, maximum=12)
    
    # change your urls to the dask dashboard so that you can see it
    hostname = socket.getfqdn()
    
    dask.config.set({'distributed.dashboard.link':'https://ondemand.hpc.ucar.edu/stable/user/{USER}/proxy/{port}/status'})
    dask.config.set({"distributed.dashboard.link":
            f"https://ondemand.hpc.ucar.edu/rnode/{hostname}/{{port}}/status"
    })
    
    # Setup your client
    client = Client(cluster)

    return cluster
