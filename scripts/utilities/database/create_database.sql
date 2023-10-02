ATTACH DATABASE file_name AS database_name;

CREATE TABLE database_name.LosAlamos_seq (
  accession_id TEXT PRIMARY KEY,
  sequence_name TEXT,
  subtype TEXT
);
