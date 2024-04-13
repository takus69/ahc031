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

    for rd in &mut rect {
        for k in 0..input.n {
            rd.push((k, 0, k+1, input.w));
        }
    }
    for rd in rect {
        for r in rd {
            println!("{} {} {} {}", r.0, r.1, r.2, r.3);
        }
    }
}
