use proconio::input;

fn main() {
    input! {
        w: u16, d: u8, n: u8,
        a: [[i32; n]; d],
    }
    println!("{} {} {}", w, d, n);
    println!("{:?}", a);
}
