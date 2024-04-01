use proconio::input;
use std::fmt;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::collections::HashSet;

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

#[derive(Debug, Clone)]
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

    fn stack(&self, i: usize, j: usize, w: usize, a: &Vec<(usize, usize)>) -> Vec<(usize, usize, i64, usize, usize)> {
        // eprintln!("stack (i, j): ({}, {}), w: {}, len a: {}, a: {:?}", i, j, w, a.len(), a);
        let mut max_h = 0;
        let mut data: Vec::<(usize, usize, i64, usize, usize)> = Vec::with_capacity(a.len());  // (h: 高さ, a: 予約の面積, b: 不足面積, k: 予約の面積のインデックス)
        for (k, ai) in a {
            let mut h = ai/w;
            if h == 0 {
                h = 1;
            }
            let b = *ai as i64 - (w * h) as i64;
            if max_h + h > self.w - a.len() {
                if self.w > a.len() + max_h {
                    h = self.w - a.len() - max_h;
                } else {
                    h = 1;
                }
            }
            data.push((h, *ai, b, *k, h));
            max_h += h;
        }

        // 不足分を調整
        data.sort_by_key(|&(_, _, b, _, _)| std::cmp::Reverse(b));  // 不足分が大きい方から処理
        let mut temp_data: Vec<(usize, usize, i64, usize, usize)> = Vec::with_capacity(data.len());
        // eprintln!("{} {}", self.w, max_h);
        let hh = 1;
        // let hh = (self.w - max_h - 1) / data.len() + 1;
        for (h, a, b, k, min_h) in &data {
            if max_h >= self.w {
                temp_data.push((*h, *a, *b, *k, *min_h));
                continue;
            }
            let mut h = h + hh;
            max_h += hh;
            if max_h > self.w {
                h -= max_h - self.w;
            }
            temp_data.push((h, *a, *a as i64 -(self.w*h) as i64, *k, *min_h));
        }

        data.copy_from_slice(&temp_data[..(temp_data.len())]);

        data
    }

    fn conv_rooms(&self, i: usize, j: usize, w: usize, data: &Vec<(usize, usize, i64, usize, usize)>) -> Vec<(Room, usize)> {
        // Roomに変換
        let mut ii = i;
        let jj = j;
        let mut rooms: Vec::<(Room, usize)> = Vec::with_capacity(data.len());  // (Room, k: 予約の面積のインデックス)
        for (h, a, _, k, min_h) in data {
            // eprintln!("top_left: ({}, {}), bottom_right: ({}, {})", ii, jj, ii+h, jj+w);
            let room = Room::new((ii, jj), (ii+h, jj+w), self.w, *a);
            // eprintln!("room: {}", room);
            rooms.push((room, *k));
            ii += h;
        }

        rooms
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

#[derive(Debug, Clone)]
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
            //cost -= pre_ar.overlap(ar) * 2;
            // eprintln!("d: {}, pre_ar overlap: {}", d, pre_ar.overlap(ar));
            pre_ar = ar;
            // eprintln!("d: {}, cost: {} {}", d, cost, ar.cost());
        }
        cost
    }

    fn optimize(&self, data1: Vec<(usize, usize, i64, usize, usize)>, data2: Vec<(usize, usize, i64, usize, usize)>) -> (Vec<(usize, usize, i64, usize, usize)>, Vec<(usize, usize, i64, usize, usize)>) {
        let mut max_h1 = 0;
        let mut max_h2 = 0;
        let mut candidates: Vec<(i64, usize, usize)> = Vec::new();
        let mut h_set: HashSet<i64> = HashSet::new();
        for (i, (h1, b1, c1, k1, min_h1)) in data1.iter().enumerate() {
            max_h1 += h1;
            for (j, (h2, b2, c2, k2, min_h2)) in data2.iter().enumerate() {
                if i == 0 {
                    max_h2 += h2;
                }
                let h: i64 = *h1 as i64 - *h2 as i64;
                candidates.push((h, i, j));
                h_set.insert(h.abs());
            }
        }
        let mut h_set2: Vec<i64> = Vec::new();
        for h in h_set {
            h_set2.push(h);
        }
        h_set2.sort_by_key(|&h| h);
        // eprintln!("{:?}", h_set2);

        let mut optimal_data1: Vec<(usize, usize, i64, usize, usize)> = Vec::new();
        let mut optimal_data2: Vec<(usize, usize, i64, usize, usize)> = Vec::new();
        let mut flag1 = vec![false; data1.len()];
        let mut flag2 = vec![false; data2.len()];
        let w: usize = 1000;
        for hhh in h_set2 {
            let hh: usize = hhh as usize;
            if max_h1 + hh > w && max_h2 + hh > w {
                break;
            }
            for (h, i, j) in &candidates {
                if h.abs() as usize == hh && !flag1[*i] && !flag2[*j] {
                    let (h1, b1, c1, k1, min_h1) = data1[*i];
                    let (h2, b2, c2, k2, min_h2) = data2[*j];
                    if *h > 0 {
                        if max_h1 + hh > w {
                            break;
                        }
                        optimal_data1.push((h1+hh, b1, c1, k1, min_h1));
                        optimal_data2.push((h2, b2, c2, k2, min_h2));
                        max_h1 += hh;
                    } else {
                        if max_h2 + hh > w {
                            break;
                        }
                        optimal_data1.push((h1, b1, c1, k1, min_h1));
                        optimal_data2.push((h2+hh, b2, c2, k2, min_h2));
                        max_h2 += hh;
                    }
                    flag1[*i] = true;
                    flag2[*j] = true;
                }
            }
        }
        for i in 0..(data1.len()) {
            if !flag1[i] {
                optimal_data1.push(data1[i].clone());
            }
        }
        for i in 0..(data2.len()) {
            if !flag2[i] {
                optimal_data2.push(data2[i].clone());
            }
        }
        (optimal_data1, optimal_data2)
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

    fn recon(&self, data1: Vec<(usize, usize, i64, usize, usize)>, data2: Vec<(usize, usize, i64, usize, usize)>, w2: usize, d: usize) -> Arrangement {
        let mut ar = self.solution.arrangements[d].clone();
        let rooms1 = ar.conv_rooms(0, 0, self.w-w2, &data1);
        let rooms2 = ar.conv_rooms(0, self.w-w2, w2, &data2);
        for (room, k) in rooms1 {
            ar.rooms[k] = room;
        }
        for (room, k) in rooms2 {
            ar.rooms[k] = room;
        }
        ar
    }

    fn solve(&mut self) {
        // 初期配置の決定
        let mut optimal_data1: Vec<Vec<(usize, usize, i64, usize, usize)>> = Vec::with_capacity(self.d);
        let mut optimal_data2: Vec<Vec<(usize, usize, i64, usize, usize)>> = Vec::with_capacity(self.d);
        let mut optimal_param: Vec<usize> = Vec::with_capacity(self.d);
        for d in 0..self.d {
            // 初期配置
            let mut ar = Arrangement::new(self.w, self.n);
            // ar.init(self.a[d].clone());
            let a: Vec<(usize, usize)> = self.a[d].iter().enumerate().map(|(i, &x)| (i, x)).collect();

            // 一列に配置
            let data = ar.stack(0, 0, self.w, &a);
            let mut rooms = ar.conv_rooms(0, 0, self.w, &data);
            rooms.sort_by_key(|&(_, k)| k);
            for (room, _) in rooms {
                ar.rooms.push(room);
            }
            let mut min_cost = ar.cost();
            let mut optimal_ar = ar.clone();
            optimal_data1.push(data.to_vec());
            optimal_data2.push(data.to_vec());
            optimal_param.push(0);

            // 二列に配置。何個か一列目から除外して、二列目に配置
            for i in 1..(self.n) {
                let sum_a: usize = a.iter().map(|(_, x)| x).take(i).sum();
                for j in 1..5 {
                    let w2 = sum_a / self.w + j;
                    if w2 > self.w / 2 { break; }
                    // 一列目
                    let data1 = ar.stack(0, 0, self.w-w2, &a[i..].to_vec());
                    let rooms1 = ar.conv_rooms(0, 0, self.w-w2, &data1);
                    // 二列目
                    let data2 = ar.stack(0, self.w-w2, w2, &a[..i].to_vec());
                    let rooms2 = ar.conv_rooms(0, self.w-w2, w2, &data2);
                    for (room, k) in rooms1 {
                        ar.rooms[k] = room;
                    }
                    for (room, k) in rooms2 {
                        ar.rooms[k] = room;
                    }
                    if ar.cost() < min_cost {
                        min_cost = ar.cost();
                        optimal_ar = ar.clone();
                        optimal_data1[d] = data1.to_vec();
                        optimal_data2[d] = data2.to_vec();
                        optimal_param[d] = w2;
                    }
                }
            }

            self.solution.arrangements.push(optimal_ar);
        }

        // 同じ高さの部屋があれば上に詰めていく
        let mut ar = Arrangement::new(self.w, self.n);
        for d in (1..self.d).rev().step_by(2) {
            let data1_1 = optimal_data1[d].clone();
            let data1_2 = optimal_data1[d-1].clone();
            let (optimal_data1_1, optimal_data1_2) = self.solution.optimize(data1_1, data1_2);
            let data2_1 = optimal_data2[d].clone();
            let data2_2 = optimal_data2[d-1].clone();
            let (optimal_data2_1, optimal_data2_2) = self.solution.optimize(data2_1, data2_2);
            optimal_data1[d] = optimal_data1_1.clone();
            optimal_data2[d] = optimal_data2_1.clone();
            optimal_data1[d-1] = optimal_data1_2.clone();
            optimal_data2[d-1] = optimal_data2_2.clone();
        
            // 部屋の再構築
            let w2_1 = optimal_param[d];
            let ar = self.recon(optimal_data1_1, optimal_data2_1, w2_1, d);
            self.solution.arrangements[d] = ar;
            let w2_2 = optimal_param[d-1];
            let ar = self.recon(optimal_data1_2, optimal_data2_2, w2_2, d-1);
            self.solution.arrangements[d-1] = ar;
        }

        let mut optimal_cost = self.solution.cost();
        for d in 0..self.d {
            let mut ar = self.solution.arrangements[d].clone();
            let mut data1 = optimal_data1[d].clone();
            let mut data2 = optimal_data2[d].clone();
            // data1.sort_by_key(|&(h, _, _, _, _)| h);
            // data2.sort_by_key(|&(h, _, _, _, _)| h);

            // 部屋をフルに使う
            let mut max_h = 0;
            for (h, _, _, _, _,) in &data1 {
                max_h += h;
            }
            let n = data1.len();
            data1[n-1].0 += self.w - max_h;
            for (h, _, _, _, _,) in &data2 {
                max_h += h;
            }
            let mut max_h = 0;
            for (h, _, _, _, _,) in &data2 {
                max_h += h;
            }
            let n = data2.len();
            data2[n-1].0 += self.w - max_h;

            // 部屋の再構築
            let w2 = optimal_param[d];
            let ar = self.recon(data1, data2, w2, d);
            self.solution.arrangements[d] = ar;
        }

        // ゼロコスト判定
        self.check_zero();
    }

    fn ans(&self) {
        for ar in self.solution.arrangements.iter() {
            for room in ar.rooms.iter() {
                println!("{} {} {} {}", room.top_left.0, room.top_left.1, room.bottom_right.0, room.bottom_right.1);
            }
        }
    }

    fn result(&self) {
        let sums: Vec<usize> = self.a.iter().map(|row| row.iter().sum()).collect();
        let average: f64 = sums.iter().sum::<usize>() as f64 / sums.len() as f64;
        let e = (1.0 - average / (self.w * self.w) as f64).sqrt();
        eprintln!("{{ \"d\": {}, \"n\": {}, \"e\": {}, \"cost\": {} }}", self.d, self.n, e, self.solution.cost());
    }

    fn check_zero(&mut self) {
        // ゼロコスト判定
        let mut max_a: Vec<usize> = vec![0; self.n];
        for i in 0..self.n {
            for d in 0..self.d {
                max_a[i] = max_a[i].max(self.a[d][i]);
            }
        }
        let mut sum_a = 0;
        for a in &max_a {
            sum_a += a;
        }
        if sum_a <= self.w*self.w {
            eprintln!("Ok zero: sum a: {}, {:?}", sum_a, max_a);
            let mut ar = Arrangement::new(self.w, self.n);
            let a: Vec<(usize, usize)> = max_a.iter().enumerate().map(|(i, &x)| (i, x)).collect();
            // 一列に配置
            let data = ar.stack(0, 0, self.w, &a);
            let mut rooms = ar.conv_rooms(0, 0, self.w, &data);
            rooms.sort_by_key(|&(_, k)| k);
            for (room, _) in rooms {
                ar.rooms.push(room);
            }
            let mut min_cost = ar.cost();
            let mut optimal_ar = ar.clone();

            // 二列に配置。何個か一列目から除外して、二列目に配置
            for i in 1..(self.n) {
                for j in 1..5 {
                    let w2 = sum_a / self.w + j;
                    if w2 > self.w / 2 { break; }
                    // 一列目
                    let data1 = ar.stack(0, 0, self.w-w2, &a[i..].to_vec());
                    let rooms1 = ar.conv_rooms(0, 0, self.w-w2, &data1);
                    // 二列目
                    let data2 = ar.stack(0, self.w-w2, w2, &a[..i].to_vec());
                    let rooms2 = ar.conv_rooms(0, self.w-w2, w2, &data2);
                    for (room, k) in rooms1 {
                        ar.rooms[k] = room;
                    }
                    for (room, k) in rooms2 {
                        ar.rooms[k] = room;
                    }
                    if ar.cost() < min_cost {
                        min_cost = ar.cost();
                        optimal_ar = ar.clone();
                    }
                }
            }
            for d in 0..self.d {
                self.solution.arrangements[d] = optimal_ar.clone();
            }
        }
    }
}

fn main() {
    let input: Input = parse_input();
    let mut solver = Solver::new(input);
    solver.solve();
    solver.ans();
    solver.result();
}
