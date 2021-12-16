use methods::{adams, f1, f2};
fn main() {
    let func_list: Vec<&dyn Fn(f64, f64, f64) -> f64> = vec![&f1, &f2];
    let answer_adams = adams(2.0, 4.0, 10, func_list.clone(), &mut vec![1.0, -1.0]);
    print_answer(answer_adams, "Адамсу");
}
fn print_answer(answer: Vec<(f64, f64, f64)>, name: &str) {
    println!("-----------------------------------------------------------------");
    println!("Ответы по {}", name);
    for i in 0..answer.len() {
        println!(
            "x: {:.02}, y1: {:.05}, y2: {:.05}",
            answer[i].0, answer[i].1, answer[i].2
        );
    }
}
