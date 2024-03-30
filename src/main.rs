use proconio::input;
use std::fmt;
use rand::seq::SliceRandom;

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

#[derive(Debug, Clone)]
struct Room {
    top_left: (usize, usize),
    bottom_right: (usize, usize),
    w: usize,
    reserved: usize,
}

impl Room {
    fn new(top_left: (usize, usize), bottom_right: (usize, usize), w: usize, reserved: usize) -> Self {
        assert!(top_left.0 <= bottom_right.0 && top_left.1 <= bottom_right.1, "Invalid Room coordinates");
        assert!(top_left.0 <= w && top_left.1 <= w, "Invalid top_left");
        assert!(bottom_right.0 <= w && bottom_right.1 <= w, "Invalid bottom_right");
        Room { top_left, bottom_right, w, reserved }
    }

    fn area(&self) -> usize {
        let width = self.top_left.0.wrapping_sub(self.bottom_right.0);
        let height = self.top_left.1.wrapping_sub(self.bottom_right.1);
        width.wrapping_mul(height)
    }

    fn area_cost(&self) -> usize {
        let mut cost = 0;
        if self.reserved > self.area() {
            cost += 100 * (self.reserved - self.area());
        }
        cost
    }

    fn partition_cost(&self) -> usize {
        let mut cost = 0;

        cost += (self.bottom_right.0 - self.top_left.0) * 2;
        cost += (self.bottom_right.1 - self.top_left.1) * 2;

        // 外周はコスト減
        if self.top_left.0 == 0 {
            cost -= self.bottom_right.1 - self.top_left.1;
        }
        if self.top_left.1 == 0 {
            cost -= self.bottom_right.0 - self.top_left.0;
        }
        if self.bottom_right.0 == self.w {
            cost -= self.bottom_right.1 - self.top_left.1;
        }
        if self.bottom_right.1 == self.w {
            cost -= self.bottom_right.0 - self.top_left.0;
        }
        cost
    }

    fn cost(&self) -> usize {
        self.area_cost() + self.partition_cost()
    }

    fn overlap_partition(&self, other: &Room) -> usize {
        let mut adjacent = 0;

        // 上下の重なり
        let mut horizontal_adjacent = self.bottom_right.1.min(other.bottom_right.1).wrapping_sub(self.top_left.1.max(other.top_left.1));
        if horizontal_adjacent > self.w {
            horizontal_adjacent = 0;
        }

        adjacent += 
            if self.bottom_right.0 == other.top_left.0 || other.bottom_right.0 == self.top_left.0 {
                horizontal_adjacent
            } else { 0 };

        // 上上、下下の重なり
        adjacent +=
            if self.top_left.0 == other.top_left.0 && self.top_left.0 != 0 {
                horizontal_adjacent
            } else { 0 };
        adjacent +=
            if self.bottom_right.0 == other.bottom_right.0 && self.bottom_right.0 != self.w {
                horizontal_adjacent
            } else { 0 };

        // 左右の重なり
        let mut vertical_adjacent = self.bottom_right.0.min(other.bottom_right.0).wrapping_sub(self.top_left.0.max(other.top_left.0));
        if vertical_adjacent > self.w {
            vertical_adjacent = 0;
        }
        adjacent += 
            if self.bottom_right.1 == other.top_left.1 || other.bottom_right.1 == self.top_left.1 {
                vertical_adjacent
            } else { 0 };

        // 左左、右右の重なり
        adjacent +=
            if self.top_left.1 == other.top_left.1 && self.top_left.1 != 0 {
                vertical_adjacent
            } else { 0 };
        adjacent +=
            if self.bottom_right.1 == other.bottom_right.1 && self.bottom_right.1 != self.w {
                vertical_adjacent
            } else { 0 };

        adjacent
    }

    fn check_adjacent(&self, other: &Room) -> (bool, usize) {
        // 二つのRoomが水平方向にかぶっているかをチェック
        let horizontal_overlap = 
            self.bottom_right.0 > other.top_left.0 && other.bottom_right.0 > self.top_left.0;

        // 二つのRoomが垂直方向にかぶっているかをチェック
        let vertical_overlap = 
            self.bottom_right.1 > other.top_left.1 && other.bottom_right.1 > self.top_left.1;

        // 二つのRoomがかぶっているかどうかを判定
        let overlap = horizontal_overlap && vertical_overlap;

        // 接している場合は、接している長さを計算
        let horizontal_adjacent = 
            if self.bottom_right.0 == other.top_left.0 || other.bottom_right.0 == self.top_left.0 {
                let mut ret = self.bottom_right.1.min(other.bottom_right.1).wrapping_sub(self.top_left.1.max(other.top_left.1));
                if ret > self.w {
                    ret = 0;
                }
                ret
            } else {
                0
            };

        let vertical_adjacent = 
            if self.bottom_right.1 == other.top_left.1 || other.bottom_right.1 == self.top_left.1 {
                let mut ret = self.bottom_right.0.min(other.bottom_right.0).wrapping_sub(self.top_left.0.max(other.top_left.0));
                if ret > self.w {
                    ret = 0;
                }
                ret
            } else {
                0
            };
        let adjacent = horizontal_adjacent + vertical_adjacent;

        (overlap, adjacent)
        
    }

    fn add_coordinate(&self, x: usize, d: i64) -> usize {
        let ret = x as i64 + d;
        if ret >= 0 && ret <= self.w as i64 {
            ret as usize
        } else {
            x
        }
    }

    fn expand(&self, d: char) -> Room {
        let mut room = self.clone();
        match d {
            'L' => room.top_left.0 = self.add_coordinate(room.top_left.0, -1),
            'R' => room.bottom_right.0 = self.add_coordinate(room.bottom_right.0, 1),
            'T' => room.top_left.1 = self.add_coordinate(room.top_left.1, -1),
            'B' => room.bottom_right.1 = self.add_coordinate(room.bottom_right.1, 1),
            _ => {},
        }
        room
    }

    fn shrink(&self, d: char) -> Room {
        let mut room = self.clone();
        match d {
            'L' => room.top_left.0 = self.add_coordinate(room.top_left.0, 1),
            'R' => room.bottom_right.0 = self.add_coordinate(room.bottom_right.0, -1),
            'T' => room.top_left.1 = self.add_coordinate(room.top_left.1, 1),
            'B' => room.bottom_right.1 = self.add_coordinate(room.bottom_right.1, -1),
            _ => {},
        }
        room
    }
    
    fn shift(&self, d: char) -> Room {
        let mut room = self.clone();
        match d {
            'L' => {
                room.top_left.0 = self.add_coordinate(room.top_left.0, -1);
                room.bottom_right.0 = self.add_coordinate(room.bottom_right.0, -1);
            },
            'R' => {
                room.top_left.0 = self.add_coordinate(room.top_left.0, 1);
                room.bottom_right.0 = self.add_coordinate(room.bottom_right.0, 1);
            },
            'T' => {
                room.top_left.1 = self.add_coordinate(room.top_left.1, -1);
                room.bottom_right.1 = self.add_coordinate(room.bottom_right.1, -1);
            },
            'B' => {
                room.top_left.1 = self.add_coordinate(room.top_left.1, 1);
                room.bottom_right.1 = self.add_coordinate(room.bottom_right.1, 1);
            },
            _ => {},
        }
        room
    }

    fn candidates(&self) -> Vec<(usize, Room)> {
        let dir = ['L', 'R', 'T', 'B'];
        let mut candidates = Vec::new();
        // 移動
        let d = dir.choose(&mut rand::thread_rng()).unwrap();
        let room = self.shift(*d);
        candidates.push((room.cost(), room));
        //for d in dir {
        //    let room = self.shift(d);
        //    candidates.push((room.cost(), room));
        //}
        // 拡張
        let d = dir.choose(&mut rand::thread_rng()).unwrap();
        let room = self.expand(*d);
        candidates.push((room.cost(), room));
        //for d in dir {
        //    let room = self.expand(d);
        //    candidates.push((room.cost(), room));
        //}
        // 縮小
        let d = dir.choose(&mut rand::thread_rng()).unwrap();
        let room = self.shrink(*d);
        candidates.push((room.cost(), room));
        //for d in dir {
        //    let room = self.shrink(d);
        //    candidates.push((room.cost(), room));
        //}
        candidates.push((self.cost(), self.clone()));
        
        // コストの昇順
        candidates.sort_by_key(|&(num, _)| num);

        candidates
    }
}

impl fmt::Display for Room {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f, "top_left: ({}, {}) bottom_right: ({}, {}), reserved: {}, area: {}",
            self.top_left.0, self.top_left.1, self.bottom_right.0, self.bottom_right.1, self.reserved, self.area())?;
        Ok(())
    }
}

#[derive(Debug)]
struct Arrangement {
    w: usize,
    n: usize,
    rooms: Vec<Room>,
}

impl Arrangement {
    fn new(w: usize, n: usize) -> Self {
        let rooms = Vec::with_capacity(n);
        Arrangement { w, n, rooms }
    }

    fn init(&mut self, a: Vec<usize>) {
        let col_n = (self.n as f32).sqrt().ceil() as usize;
        let room_len = self.w / col_n;
        for (i, ai) in a.iter().enumerate().take(self.n) {
            let li = (i / col_n)*room_len;
            let lj = (i % col_n)*room_len;
            let room = Room::new((li, lj), (li+room_len, lj+room_len), self.w, *ai);
            self.rooms.push(room);
        }
    }

    fn optimize(&mut self) {
        for i in 0..self.n {
            let candidates = self.rooms[i].candidates();
            let mut min_cost= std::usize::MAX;
            let mut optimal_candidate = None;
            // println!("optimize {}-room: {}, cost: {}", i, self.rooms[i], self.rooms[i].cost());

            for (_, candidate) in &candidates {
                // println!("candidate: {:?}, cost: {}", candidate, candidate.cost());
                let mut candidate_overlap = false;

                for j in 0..self.n {
                    if i == j { continue; }
                    let other = &self.rooms[j];
                    let (overlap, _) = candidate.check_adjacent(other);
                    // println!("overlap check {} adjacent: {} other: {}", overlap, adjacent, other);

                    if overlap {
                        candidate_overlap = true;
                        // println!("overlap! {}", other);
                        break; // 重なっている場合は次のcandidateへ
                    }
                }

                if !candidate_overlap {
                    // println!("final candidate: {}, cost: {}, adjacent: {}", candidate, candidate.cost(), candidate_adjacent);
                    let cost = candidate.cost();
                    if cost < min_cost {
                        min_cost = cost;
                        optimal_candidate = Some(candidate);
                    }
                }
            }
            // println!("candidate: {:?}", optimal_candidate);
            if let Some(room) = optimal_candidate {
                // println!("room modify: {}", room);
                self.rooms[i] = room.clone();
            }
        }
    }

    fn area_cost(&self) -> usize {
        let mut cost = 0;
        for i in 0..self.n {
            let room1 = &self.rooms[i];
            cost += room1.area_cost();
        }
        cost
    }

    fn partition_cost(&self) -> usize {
        let mut cost = 0;
        for i in 0..self.n {
            let room1 = &self.rooms[i];
            cost += room1.partition_cost();
            for j in (i+1)..self.n {
                let room2 = &self.rooms[j];
                let (overlap, adjacent) = room1.check_adjacent(room2);
                if overlap {
                    eprintln!("overlap! room1: {}, room2: {}", room1, room2);
                }
                cost -= adjacent;
            }
        }
        cost
    }

    fn cost(&self) -> usize {
        self.area_cost() + self.partition_cost()
    }

    fn overlap(&self, ar: &Arrangement) -> usize {
        let mut adjacent = 0;
        for room1 in &self.rooms {
            for room2 in &ar.rooms {
                adjacent += room1.overlap_partition(room2);
            }
        }
        adjacent
    }
}

impl fmt::Display for Arrangement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for room in &self.rooms {
            writeln!(f, "{}", room)?;
        }
        Ok(())
    }
}

#[derive(Debug)]
struct Solution {
    d: usize,
    arrangements: Vec<Arrangement>,
}

impl Solution {
    fn new(d: usize) -> Self {
        let arrangements = Vec::with_capacity(d);
        Solution { d, arrangements }
    }

    fn cost(&self) -> usize {
        let mut cost = 1;
        let mut pre_ar = &self.arrangements[0];

        cost += pre_ar.area_cost();
        // eprintln!("d: 0, cost: {}", cost);
        for d in 1..self.d {
            let ar = &self.arrangements[d];
            cost += ar.cost();
            // eprintln!("d: {}, ar cost: {}", d, ar.cost());
            // パーティションの撤去
            cost += pre_ar.partition_cost();
            // eprintln!("d: {}, partition cost: {}", d, pre_ar.partition_cost());
            // 共通部分はコスト減(撤去と設置分)
            cost -= pre_ar.overlap(ar) * 2;
            // eprintln!("d: {}, pre_ar overlap: {}", d, pre_ar.overlap(ar));
            pre_ar = ar;
            // eprintln!("d: {}, cost: {} {}", d, cost, ar.cost());
        }
        cost
    }

}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for ar in &self.arrangements {
            writeln!(f, "{}", ar)?;
        }
        Ok(())
    }
}

struct Solver {
    w: usize,
    d: usize,
    n: usize,
    a: Vec<Vec<usize>>,
    solution: Solution,
}

impl Solver {
    fn new(input: Input) -> Solver {
        Solver {
            w: input.w,
            d: input.d,
            n: input.n,
            a: input.a,
            solution: Solution::new(input.d),
        }
    }

    fn solve(&mut self) {
        // 初期配置の決定
        for d in 0..self.d {
            // 初期化
            let mut ar = Arrangement::new(self.w, self.n);
            ar.init(self.a[d].clone());
            self.solution.arrangements.push(ar);
        }

        // 最適化
        for ar in &mut self.solution.arrangements[..self.d] {
            for _ in 0..3000 {
                ar.optimize();
                // println!("{}", ar);
            }
        }

        // 日付間の配置調整
    }

    fn ans(&self) {
        for ar in self.solution.arrangements.iter() {
            for room in ar.rooms.iter() {
                println!("{} {} {} {}", room.top_left.0, room.top_left.1, room.bottom_right.0, room.bottom_right.1);
            }
        }
    eprintln!("cost: {}", self.solution.cost());
    }
}

fn main() {
    let input: Input = parse_input();
    let mut solver = Solver::new(input);
    solver.solve();
    solver.ans();
}
