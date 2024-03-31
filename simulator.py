import numpy as np
from tqdm import tqdm
import pandas as pd
import time
import multiprocessing
import subprocess
import json

def run(i):
    output_str = subprocess.run(f'powershell cat in/{i:04}.txt | .\\target\\debug\\ahc031.exe > out/{i:04}.txt', shell=True, capture_output=True, text=True).stderr
    result = json.loads(output_str.split('\n')[0])
    return result

def main(i):
    start = time.time()
    # print(i, 'start')
    r = run(i)
    t = round(time.time()-start, 4)
    d = r['d']
    n = r['n']
    e = r['e']
    cost = r['cost']
    cost2 = r['cost'] / ((d-1)*(n-2))
    data = [i, d, n, round(e, 4), cost, cost2, t]
    print('\r', 'end', i, end='')
    # print(i, 'end')
    return data


if __name__ == '__main__':
    start = time.time()
    trial = 200
    result = []
    '''
    for i in tqdm(range(trial)):
        r = run(i)
        t = round(time.time()-start, 4)
        d = r['d']
        n = r['n']
        e = r['e']
        cost = r['cost']
        data = [i, d, n, round(e, 4), cost, t]
        result.append(data)
    '''
    processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=processes) as pool:
        data = [pool.apply_async(main, (i,)) for i in range(trial)]
        result = [d.get() for d in data]
    print()
    df = pd.DataFrame(result, columns=['i', 'd', 'n', 'e', 'cost', 'cost2', 'time'])
    cost = np.mean(df['cost']) * 50
    cost2 = np.mean(df['cost2']) * 50
    print(f'cost:', format(int(df['cost'].sum()*50/trial), ','), 'cost2:', format(int(df['cost2'].mean()), ','))
    df.to_csv('result.csv', index=False)
    print(f'end elapsed time: {time.time()-start:.2f}s')
