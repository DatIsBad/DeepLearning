CREATE TABLE IF NOT EXISTS samples (
    id INTEGER PRIMARY KEY,
    enzyme TEXT NOT NULL,
    sample TEXT,
    orf TEXT,
    rec_sequence TEXT,
    size INTEGER NOT NULL,
    fragment int NOT NULL,
    filename TEXT NOT NULL,
    line INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS alignments (
    id INTEGER PRIMARY KEY AUTOINCREMENT,

    enzyme_a_id INTEGER NOT NULL,
    enzyme_b_id INTEGER NOT NULL,

    algorithm TEXT NOT NULL,
    score INTEGER NOT NULL,

    aligned_seq_a TEXT NOT NULL,
    aligned_seq_b TEXT NOT NULL,

    UNIQUE(enzyme_a_id, enzyme_b_id, algorithm) ON CONFLICT REPLACE,

    FOREIGN KEY (enzyme_a_id) REFERENCES samples(id),
    FOREIGN KEY (enzyme_b_id) REFERENCES samples(id)
);
