use proconio::input;
use std::fmt;

struct Input {
    w: u16,
    d: u16,
    n: u16,
    a: Vec<Vec<u32>>,
}

fn parse_input() -> Input {
    input! {
        w: u16, d: u16, n: u16,
        a: [[u32; n]; d],
    }
    Input { w, d, n, a }
}

struct Room {
    top_left: (u16, u16),
    bottom_right: (u16, u16),
}

impl Room {
    fn new(top_left: (u16, u16), bottom_right: (u16, u16)) -> Room {
        assert!(top_left.0 <= bottom_right.0 && top_left.1 <= bottom_right.1,
        "Invalid Room coordinates");

        Room { top_left, bottom_right }
    }

    fn area(&self) -> u16 {
        (self.top_left.0 - self.bottom_right.0)*(self.top_left.1 - self.bottom_right.1)
    }
}

impl fmt::Display for Room {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {} {}", self.top_left.0, self.top_left.1, self.bottom_right.0, self.bottom_right.1)?;
        Ok(())
    }
}

#[derive(Default)]
struct Arrangement {
    rooms: Vec<Room>,
}

impl fmt::Display for Arrangement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for room in self.rooms.iter() {
            writeln!(f, "{}", room)?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct Solution {
    arrangements: Vec<Arrangement>,
}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for ar in self.arrangements.iter() {
            write!(f, "{}", ar)?;
        }
        Ok(())
    }
}

struct Solver {
    w: u16,
    d: u16,
    n: u16,
    a: Vec<Vec<u32>>,
    solution: Solution,
}

impl Solver {
    fn new(input: Input) -> Solver {
        Solver {
            w: input.w,
            d: input.d,
            n: input.n,
            a: input.a,
            solution: Solution::default(),
        }
    }

    fn solve(&mut self) {
        for _ in 0..self.d {
            let mut ar = Arrangement::default();
            for j in 0..self.n {
                ar.rooms.push(Room::new((j, 0), (j+1, self.w)));
            }
            self.solution.arrangements.push(ar);
        }

    }

    fn ans(&self) {
        println!("{}", self.solution);
    }
}

fn main() {
    let input: Input = parse_input();
    let mut solver = Solver::new(input);
    solver.solve();
    solver.ans();
}
