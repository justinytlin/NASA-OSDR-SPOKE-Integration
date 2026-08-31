"""Python client for the SPOKE REST API (https://spoke.rbvi.ucsf.edu/swagger/).

The API is public and read-only. All endpoints return JSON.
Graph-returning endpoints (metagraph, neighborhood, expand) use a
cytoscape.js-style format: {"data": [...nodes...], "edges"/[...]}.
"""

from urllib.parse import quote

import requests

BASE_URL = "https://spoke.rbvi.ucsf.edu/api/v1"


class SpokeClient:
    def __init__(self, base_url: str = BASE_URL, timeout: int = 120):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers["Accept"] = "application/json"

    def _get(self, path: str, params: dict | None = None):
        url = f"{self.base_url}/{path}"
        resp = self.session.get(url, params=params or {}, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json()

    # ── Schema / metadata ────────────────────────────────────────────────

    def version(self) -> dict:
        """SPOKE version plus per-source database update timestamps."""
        return self._get("version")

    def metagraph(self) -> dict:
        """The SPOKE metagraph: node types and how they interconnect."""
        return self._get("metagraph")

    def types(self) -> dict:
        """Node types, edge types, and default query parameters."""
        return self._get("types")

    # ── Lookup / search ──────────────────────────────────────────────────

    def search(self, query: str, node_type: str | None = None) -> list:
        """Lucene search for nodes, optionally restricted to one node type.

        Returns a list of {node_type, identifier, name, score}.
        """
        if node_type:
            path = f"search/{quote(node_type, safe='')}/{quote(query, safe='')}"
        else:
            path = f"search/{quote(query, safe='')}"
        return self._get(path)

    def node(self, node_type: str, attribute: str, value: str) -> dict:
        """Fetch a single node by attribute match, e.g. node('Disease', 'identifier', 'DOID:2377')."""
        path = (
            f"node/{quote(node_type, safe='')}/"
            f"{quote(attribute, safe='')}/{quote(str(value), safe='')}"
        )
        return self._get(path)

    # ── Graph traversal ──────────────────────────────────────────────────

    def neighborhood(
        self,
        node_type: str,
        attribute: str,
        value: str,
        depth: int | None = None,
        node_filters: list[str] | None = None,
        edge_filters: list[str] | None = None,
        **cutoffs,
    ) -> dict:
        """Neighborhood graph around a node.

        node_filters / edge_filters restrict which node and edge types are
        returned. Extra keyword args are passed through as cutoff params,
        e.g. cutoff_DaG_textmining=3.0, cutoff_CtD_phase=3.
        """
        params: dict = dict(cutoffs)
        if depth is not None:
            params["depth"] = depth
        if node_filters:
            params["node_filters"] = node_filters
        if edge_filters:
            params["edge_filters"] = edge_filters
        path = (
            f"neighborhood/{quote(node_type, safe='')}/"
            f"{quote(attribute, safe='')}/{quote(str(value), safe='')}"
        )
        return self._get(path, params)

    def expand(self, node_type: str, node_id: int, node_ids: list[int] | None = None, **params) -> dict:
        """Expand a node by its internal SPOKE id (from a previous graph result)."""
        if node_ids:
            params["node_ids"] = node_ids
        return self._get(f"expand/{quote(node_type, safe='')}/{node_id}", params)

    def sea(self, smiles_or_zinc: str) -> dict:
        """Similarity Ensemble Approach search by SMILES string or ZINC id."""
        return self._get(f"sea/{quote(smiles_or_zinc, safe='')}")


# ── Helpers for working with graph responses ─────────────────────────────

def graph_elements(graph: dict) -> list[dict]:
    """Flatten a cytoscape.js-style graph response into element dicts."""
    if isinstance(graph, list):
        return graph
    return graph.get("elements", graph.get("data", []))


def split_graph(graph: dict) -> tuple[list[dict], list[dict]]:
    """Split a graph response into (nodes, edges) using their data payloads."""
    nodes, edges = [], []
    for el in graph_elements(graph):
        data = el.get("data", el)
        (edges if "source" in data else nodes).append(data)
    return nodes, edges
