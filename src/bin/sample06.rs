use proconio::input;

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

struct Output {
    rect: Vec<Vec<(usize, usize, usize, usize)>>,
}

fn compute_score(input: &Input, out: &Output) -> (i64, String) {
    let (mut score, err) = compute_score_details(input, &out.rect);
    if !err.is_empty() {
        score = 0;
    }
    (score, err)
}

fn compute_score_details(
    input: &Input,
    out: &[Vec<(usize, usize, usize, usize)>],
) -> (i64, String) {
    let mut score = 0;
    let mut hs: Vec<(usize, usize, usize)> = Vec::new();
    let mut vs: Vec<(usize, usize, usize)> = Vec::new();
    for (d, _) in out.iter().enumerate() {
        for p in 0..input.n {
            for q in 0..p {
                if out[d][p].2.min(out[d][q].2) > out[d][p].0.max(out[d][q].0)
                    && out[d][p].3.min(out[d][q].3) > out[d][p].1.max(out[d][q].1)
                {
                    return (0, format!("Rectangles {} and {} overlap on day {}.", q, p, d));
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
            vs = vs2;
            continue;
        }
        // 横の撤去
        for a in &hs {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
            for b in &hs2 {
                let r = (a.2 as i64).min(b.2 as i64);
                let l = (a.1 as i64).max(b.1 as i64);
                if a.0 == b.0 && l < r{
                    score -= r - l;
                }
            }
        }
        // 縦の撤去
        for a in &vs {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
            for b in &vs2 {
                let r = (a.2 as i64).min(b.2 as i64);
                let l = (a.1 as i64).max(b.1 as i64);
                if a.0 == b.0 && l > r{
                    score -= l - r;
                }
            }
        }
        // 横の設置
        for a in &hs2 {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
        }
        // 縦の設置
        for a in &vs2 {
            if a.0 == 0 || a.0 == input.w { continue; }
            score += a.2 as i64 - a.1 as i64;
        }
        hs = hs2;
        vs = vs2;
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

fn main() {
    let input: Input = parse_input();
    let rect: Vec<Vec<(usize, usize, usize, usize)>> = vec![Vec::with_capacity(input.n); input.d];
    let mut output = Output{ rect };

    for d in 0..input.d {
        let rd = &mut output.rect[d];
        let mut k1 = 0;
        let mut k2 = 0;
        for k in 0..input.n {
            k2 += (input.a[d][k] / input.w) + 1;
            rd.push((k1, 0, k2.min(input.w), input.w));
            k1 = k2;
        }
    }
    ans(&output);
    let (cost, _) = compute_score(&input, &output);
    eprintln!("cost: {}", cost);
}
