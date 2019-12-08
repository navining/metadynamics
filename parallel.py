import multiprocessing

def update(x):
    print(x)
    print(multiprocessing.current_process())
    return 1,2

def parallel():
    num_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(num_cores)
    list=[(1,2),(3,4)]
    returns = pool.map(update, list)
    print(returns)

if __name__ == '__main__':
    parallel()