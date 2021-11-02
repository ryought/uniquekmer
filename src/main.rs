use bio::data_structures::suffix_array::{lcp, suffix_array};
use bio::io::fasta;

fn main() {
    let filename = std::env::args().nth(1).expect("no filename given");
    let reader = fasta::Reader::from_file(filename).unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        // println!("{:?}", record.seq());

        let mut text = record.seq().to_vec();
        text.push(b'$');

        let length = text.len();
        let pos = suffix_array(&text);
        let lcp = lcp(&text, &pos);
        let uncompressed = lcp.decompress();
        // length of the suffix
        let lengths: Vec<usize> = pos.iter().map(|p| length - p).collect();

        let n = pos.len();
        let uniquek: Vec<isize> = (0..n)
            .map(|i| {
                let u = if i + 1 < n {
                    std::cmp::max(uncompressed[i], uncompressed[i + 1])
                } else {
                    uncompressed[i]
                };
                let l = lengths[i] as isize;
                if u + 1 < l {
                    u + 1
                } else {
                    -1
                }
            })
            .collect();

        /*
        // debug
        for i in 0..pos.len() {
            let x = pos[i];
            println!(
                "{} {} {} {} suffix={:?}",
                pos[i],
                uncompressed[i],
                lengths[i],
                uniquek[i],
                std::str::from_utf8(&text[x..]).unwrap(),
            );
        }
        */

        let mut hist: Vec<usize> = vec![0; n];
        let mut max_u: usize = 0;
        for i in 0..n {
            let u = uniquek[i];
            if u != -1 {
                max_u = std::cmp::max(u as usize, max_u);
                hist[u as usize] += 1;
            }
        }

        for i in 1..=max_u {
            println!("{}\t{}\t{}", record.id(), i, hist[i]);
        }
    }
}
