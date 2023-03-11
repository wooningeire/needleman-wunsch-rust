mod nw;

fn main() {
    nw::needleman_wunsch(
        "a fluffy amphithere",
        "that faffy amphitheatre",
        '▓',
    );
}
