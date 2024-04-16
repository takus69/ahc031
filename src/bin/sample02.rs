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

fn main() {
    let input: Input = parse_input();
    let mut rect: Vec<Vec<(usize, usize, usize, usize)>> = vec![Vec::with_capacity(input.n); input.d];

    for d in 0..input.d {
        let rd = &mut rect[d];
        let mut k1 = 0;
        let mut k2 = 0;
        for k in 0..input.n {
            k2 += (input.a[d][k] / input.w) + 1;
            rd.push((k1, 0, k2.min(input.w), input.w));
            k1 = k2;
        }
    }
    for rd in rect {
        for r in rd {
            println!("{} {} {} {}", r.0, r.1, r.2, r.3);
        }
    }
}
