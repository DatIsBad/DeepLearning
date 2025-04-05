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

