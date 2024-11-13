import configparser
import logging
import os.path
import sqlite3
import warnings

IS_IDENTITY = 1
IS_SIMILARITY = 0

log = logging.getLogger("base")


class ResolverCache:
    def __init__(self, configuration: configparser.ConfigParser) -> None:
        """
        Initializes new or loads existing resolver cache database
        :param configuration: for determining cache directory and cache file
        """

        cache_directory: str = configuration.get("cache", "directory")
        cache_file: str = configuration.get("cache", "file")

        self.cache_path: str = os.path.join(cache_directory, cache_file)
        if not os.path.exists(self.cache_path):
            log.warning(
                f"resolver cache not found at '{self.cache_path}'. Initializing empty"
            )
            self._create_empty()

        self.conn: sqlite3.Connection = sqlite3.connect(self.cache_path)

    def _create_empty(self) -> None:
        """
        Creates an empty cache database file
        :return: None
        """
        directory_path: str = os.path.dirname(self.cache_path)
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
        with sqlite3.connect(self.cache_path) as conn:
            conn.execute(
                """CREATE TABLE "CACHE" (
                            "RESOLVER_NAME"	TEXT NOT NULL,
                            "SIGNATURE"	TEXT NOT NULL,
                            "IDENTIFIER"	TEXT NOT NULL,
                            "IS_IDENTITY"	INTEGER NOT NULL,
	                        "SIMILARITY"	REAL,
                            PRIMARY KEY("IDENTIFIER","SIGNATURE","RESOLVER_NAME"));
                        """
            )

    def fetch_identity_identifiers(
        self, resolver_name: str, signature: str
    ) -> list[str]:
        """
        Finds all identifiers already resolved on identity level and cached
        :param resolver_name: of resolver trying to identify the signature
        :param signature: of the molecule
        :return: list of all identifiers associated with the signature on identity level
        """
        with self.conn as conn:
            records = conn.execute(
                "SELECT IDENTIFIER FROM CACHE WHERE RESOLVER_NAME=? AND SIGNATURE=? AND IS_IDENTITY=?",
                (resolver_name, signature, IS_IDENTITY),
            ).fetchall()
            return [record[0] for record in records]

    def fetch_similarity_identifiers(
        self, resolver_name: str, signature: str
    ) -> list[tuple[str, float]]:
        """
        Finds all similar identifiers already resolved on similarity level and cached
        :param resolver_name: of resolver trying to identify the signature
        :param signature: of the molecule
        :return: list of all identifiers associated with the fingerprint on similarity level
        """
        with self.conn as conn:
            records = conn.execute(
                "SELECT IDENTIFIER, SIMILARITY FROM CACHE "
                "WHERE RESOLVER_NAME=? AND SIGNATURE=? AND IS_IDENTITY=?",
                (resolver_name, signature, IS_SIMILARITY),
            ).fetchall()
            return records

    def save_identity_identifiers(
        self, resolver_name: str, signature: str, identifiers: list[str]
    ) -> None:
        """
        Saves all identifiers found into the cache as identical resolutions of the signature
        :param resolver_name: of the resolver that found the identifiers
        :param signature: of which identifiers were found
        :param identifiers: to be saved
        :return: None
        """

        with self.conn as conn:
            identifiers_already_cached: set[str] = set(
                self.fetch_identity_identifiers(resolver_name, signature)
            )
            for ident in identifiers:
                if ident in identifiers_already_cached:
                    continue
                conn.execute(
                    "INSERT INTO CACHE VALUES (?, ?, ?, ?, ?)",
                    (resolver_name, signature, ident, IS_IDENTITY, None),
                )
            conn.commit()

    def save_similarity_identifiers(
        self, resolver_name: str, signature: str, identifiers: list[tuple[str, float]]
    ) -> None:
        """
        Saves all identifiers found into the cache as similar resolutions of the signature
        :param resolver_name: of the resolver that found the identifiers
        :param signature: of which identifiers were found
        :param identifiers: to be saved
        :return: None
        """

        with self.conn as conn:
            identifiers_already_cached: set[str] = set(
                _ident
                for _ident, _ in self.fetch_similarity_identifiers(
                    resolver_name, signature
                )
            )
            for _ident, _similarity in identifiers:
                if _ident in identifiers_already_cached:
                    continue
                conn.execute(
                    "INSERT INTO CACHE VALUES (?, ?, ?, ?, ?)",
                    (resolver_name, signature, _ident, IS_SIMILARITY, _similarity),
                )
            conn.commit()
