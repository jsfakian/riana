#!/usr/bin/env python3
"""
Export user usage from SQLite into CSV with columns:
    id, username, first_name, last_name, usage_minutes
Tables:
  - auth_user (Django default)
  - nffa_usage: usage_minutes per user_id
"""
import sqlite3
import csv
import argparse
import sys

def export_usage(db_path: str, csv_path: str):
    # Connect to SQLite database
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.cursor()
        # Perform join between auth_user and nffa_usage
        query = """
        SELECT u.id, u.username, u.first_name, u.last_name, n.usage_minutes
        FROM nffa_usage n
        JOIN auth_user u ON n.user_id = u.id
        ORDER BY u.id
        """
        cursor.execute(query)
        rows = cursor.fetchall()

        # Write to CSV
        with open(csv_path, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Header
            writer.writerow(["id", "username", "first_name", "last_name", "usage_minutes"])
            # Data rows
            for row in rows:
                writer.writerow(row)

        print(f"Exported {len(rows)} records to '{csv_path}'.")
    except sqlite3.Error as e:
        print(f"SQLite error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export user usage data from SQLite to CSV"
    )
    parser.add_argument(
        "db_path",
        help="Path to the SQLite database file (e.g., db.sqlite3)"
    )
    parser.add_argument(
        "-o", "--output",
        default="usage.csv",
        help="Path for output CSV file (default: usage.csv)"
    )
    args = parser.parse_args()
    export_usage(args.db_path, args.output)

