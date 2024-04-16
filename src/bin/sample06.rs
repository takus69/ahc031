use proconio::input;
use std::collections::HashSet;
use rand::prelude::*;
use std::time::{Duration, Instant};

struct Input {
    w: usize,
    d: usize,
    n: usize,
    a: Vec<Vec<usize>>,
}

fn parse_input() -> Input {
    input! {
        w: usize, d: usize, n: usize,
        a: [[usize; n]; d],
    }
    Input { w, d, n, a }
}

#[derive(Clone)]
struct Output {
    rect: Vec<Vec<(usize, usize, usize, usize)>>,
}

fn compute_score(input: &Input, out: &Output) -> (i64, String) {
    let (score, err) = compute_score_details(input, &out.rect);
    (score, err)
}

fn compute_score_details(
    input: &Input,
    out: &[Vec<(usize, usize, usize, usize)>],
) -> (i64, String) {
    let mut score = 0;
    let max_score: i64 = 1_000_000_000_000;
    let mut hs: Vec<(usize, usize, usize)> = Vec::new();
    // let mut vs: Vec<(usize, usize, usize)> = Vec::new();
    for (d, _) in out.iter().enumerate() {
        for p in 0..input.n {
            if out[d][p].2 > input.w || out[d][p].3 > input.w {
                return (max_score, format!("Over {}'s w on day {}.", p, d));
            }
            for q in 0..p {
                if out[d][p].2.min(out[d][q].2) > out[d][p].0.max(out[d][q].0)
                    && out[d][p].3.min(out[d][q].3) > out[d][p].1.max(out[d][q].1)
                {
                    return (max_score, format!("Rectangles {} and {} overlap on day {}.", q, p, d));
                }
            }
        }

        let mut hs2: Vec<(usize, usize, usize)> = Vec::new();
        let mut vs2: Vec<(usize, usize, usize)> = Vec::new();
        for k in 0..input.n {
            let (i0, j0, i1, j1) = out[d][k];
            let b = (i1 - i0) * (j1 - j0);
            if input.a[d][k] > b {
                score += 100 * (input.a[d][k] - b) as i64;
            }

            hs2.push((i0, j0, j1));
            hs2.push((i1, j0, j1));
            vs2.push((j0, i0, i1));
            vs2.push((j1, i0, i1));
        }

        if d == 0 {
            hs = hs2;
            // vs = vs2;
            continue;
        }
        // 横の撤去
        for a in &hs {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
        }
        // 縦の撤去
        /*
        for a in &vs {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
        }*/
        // 横の設置
        for a in &hs2 {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
        }
        // 縦の設置
        /*
        for a in &vs2 {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
        }*/
        hs = hs2;
        // vs = vs2;
    }
    (score + 1, String::new())
}

fn ans(output: &Output) {
    for rd in &output.rect {
        for r in rd {
            println!("{} {} {} {}", r.0, r.1, r.2, r.3);
        }
    }
}

fn decide_partition(a: &[usize], n: usize, w: usize, k: usize, vs: &[usize]) -> Vec<(usize, usize, usize, usize)> {
    let mut rect: Vec<(usize, usize, usize, usize)> = Vec::with_capacity(n);
    let mut hs = vec![0; k];
    let mut free = Vec::with_capacity(k);
    let mut v1 = &vs[0];
    for v2 in &vs[1..] {
        free.push(w*(v2-v1));
        v1 = v2;
    }

    for ai in a.iter().rev() {
        let mut target_i = 0;
        let mut min_free = w * w;
        let mut v1 = vs[0];
        for i in 0..k {
            let v2 = vs[i+1];
            let h = ai / (v2 - v1) + 1;
            if h + hs[i] > w { continue; }
            let bi = (v2 - v1) * h;
            let f = free[i];
            // eprintln!("min free: {}, f: {}, bi: {}", min_free, f, bi);
            if f >= bi {
                if min_free > f - bi {
                    min_free = f - bi;
                    target_i = i;
                }
            } else { continue; }
            v1 = v2;
        }

        let v1 = vs[target_i];
        let v2 = vs[target_i+1];
        let h = ai / (v2 - v1) + 1;
        free[target_i] = min_free;
        let h1 = hs[target_i];
        let mut h2 = h1 + h;
        if h2 > w { h2 = h1 + 1; }
        hs[target_i] = h2;
        rect.push((h1, v1, h2, v2));
    }
    rect.reverse();
    rect
}

fn rnd_split(input: &Input, rng: &mut dyn RngCore) -> (usize, Vec<usize>) {
    let k = rng.gen_range(2i32..input.n as i32) as usize;
    let mut vs = HashSet::new();
    for _ in 0..k {
        let v = rng.gen_range(1i32..=input.w as i32) as usize;
        vs.insert(v);
    }
    vs.insert(0);
    vs.insert(input.w);
    let mut vs: Vec<usize> = vs.into_iter().collect();
    vs.sort();
    let k = vs.len()-1;
    // eprintln!("rnd k: {}, vs: {:?}", k, vs);

    (k, vs)
}

fn result(input: &Input, cost: i64) {
    let sums: Vec<usize> = input.a.iter().map(|row| row.iter().sum()).collect();
    let average: f64 = sums.iter().sum::<usize>() as f64 / sums.len() as f64;
    let e = (1.0 - average / (input.w * input.w) as f64).sqrt();
    eprintln!("{{ \"d\": {}, \"n\": {}, \"e\": {}, \"cost\": {} }}", input.d, input.n, e, cost);
}

fn main() {
    let start = Instant::now();

    let input: Input = parse_input();
    let rect: Vec<Vec<(usize, usize, usize, usize)>> = vec![Vec::with_capacity(input.n); input.d];
    let mut output = Output{ rect };

    let seed = 0;
    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(seed);

    let mut best_output = output.clone();
    let mut best_cost = i64::MAX;
    for _ in 0..100000 {
        if start.elapsed() >= Duration::from_secs_f64(2.8) {
            // eprintln!("Time up!!!");
            break;
        }
        let (k, vs) = rnd_split(&input, &mut rng);
        for d in 0..input.d {
            output.rect[d] = decide_partition(&input.a[d], input.n, input.w, k, &vs);
        }
        let (cost, _) = compute_score(&input, &output);
        if best_cost > cost {
            best_cost = cost;
            best_output = output.clone();
        }
    }
    ans(&best_output);
    let (cost, _) = compute_score(&input, &best_output);
    result(&input, cost);
    // eprintln!("cost: {}, s: {}", cost, s);
}
