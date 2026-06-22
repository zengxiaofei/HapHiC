#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created at: 2026-05-22 09:43

import argparse
import gzip
import logging
import os
import re
import statistics
from argparse import Namespace
from typing import TextIO, cast

import networkx as nx

_logger = logging.getLogger(__name__)


def parse_args() -> Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed arguments as ``argparse.Namespace``.
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description='Evaluate assembly quality by comparing GFA and AGP files of two genomes.'
    )
    parser.add_argument(
        '--yz_gfa',
        required=True,
        help='Path to YZ081609 GFA file (.gfa or .gfa.gz)',
    )
    parser.add_argument(
        '--yz_agp',
        required=True,
        help='Path to YZ081609 AGP file (.agp or .agp.gz)',
    )
    parser.add_argument(
        '--poj_gfa',
        required=True,
        help='Path to POJ GFA file (.gfa or .gfa.gz)',
    )
    parser.add_argument(
        '--poj_agp',
        required=True,
        help='Path to POJ AGP file (.agp or .agp.gz)',
    )
    return parser.parse_args()


def check_file_exists(*filenames: str) -> None:
    """Check if the files exist.

    Args:
        filenames: One or more file paths to check.

    Raises:
        FileNotFoundError: If any of the files does not exist.
    """
    missing_files: list[str] = []
    for filename in filenames:
        if filename and not os.path.exists(filename):
            missing_files.append(filename)
    if missing_files:
        missing_msg = ', '.join(repr(f) for f in missing_files)
        _logger.error('The following files do not exist: %s', missing_msg)
        raise FileNotFoundError(
            f'Error: The following files do not exist: {missing_msg}'
        )


def gopen(filename: str, mode: str = 'rt') -> TextIO:
    """Open a file, whether compressed (``.gz``) or not.

    Args:
        filename: Path to the file to open. If it ends with ``.gz``,
            the file will be opened as gzip compressed.
        mode: Open mode. Only ``'rt'`` and ``'wt'`` are supported.

    Returns:
        A file object opened in text mode.

    Raises:
        ValueError: If ``mode`` is not ``'rt'`` or ``'wt'``.
    """
    if mode not in ('rt', 'wt'):
        raise ValueError(
            f"Invalid mode {mode!r}, only 'rt' and 'wt' are supported"
        )
    if filename.endswith('.gz'):
        return cast(TextIO, gzip.open(filename, mode=mode, encoding='utf-8'))
    return cast(TextIO, open(filename, mode=mode, encoding='utf-8'))


def validate_input_files(args: Namespace) -> None:
    """Validate input file extensions and check for duplicate basenames.

    Args:
        args: Parsed command-line arguments.

    Raises:
        ValueError: If any file has an invalid extension or if basenames are duplicated.
    """
    ext_map: dict[str, tuple[str, ...]] = {
        'yz_gfa': ('.gfa', '.gfa.gz'),
        'poj_gfa': ('.gfa', '.gfa.gz'),
        'yz_agp': ('.agp', '.agp.gz'),
        'poj_agp': ('.agp', '.agp.gz'),
    }
    for attr_name, valid_exts in ext_map.items():
        filepath = getattr(args, attr_name)
        if not any(filepath.endswith(ext) for ext in valid_exts):
            valid_str = ' or '.join(repr(e) for e in valid_exts)
            _logger.error(
                'Invalid file extension for --%s: %s', attr_name, repr(filepath)
            )
            raise ValueError(
                f'Error: {repr(filepath)} does not have a valid extension, '
                f'expected {valid_str}'
            )

    attr_names: list[str] = ['yz_gfa', 'yz_agp', 'poj_gfa', 'poj_agp']
    basenames: list[str] = [os.path.basename(getattr(args, a)) for a in attr_names]
    seen: dict[str, str] = {}
    duplicates: list[str] = []
    for attr, name in zip(attr_names, basenames):
        if name in seen:
            duplicates.append(f'{repr(name)} (--{seen[name]} and --{attr})')
        else:
            seen[name] = attr
    if duplicates:
        duplicates_msg = ', '.join(duplicates)
        _logger.error('Duplicate input file basenames detected: %s', duplicates_msg)
        raise ValueError(
            f'Error: Duplicate input file basenames detected: {duplicates_msg}'
        )


def parse_gfa(
    gfa_file: str,
    length_tag: str = 'LN:i:',
) -> nx.DiGraph:
    """Parse a GFA file and build an assembly graph.

    Parses ``S`` (segment) and ``L`` (link) lines from a GFA file to
    construct a directed graph.  Each segment becomes a pair of forward
    and reverse nodes, and each link becomes a directed edge with
    ``match_bases`` and ``distance`` attributes.  The distance is
    computed as ``node_length - match_bases``, clamped to zero when
    ``match_bases`` exceeds the node length.

    Args:
        gfa_file: Path to the input GFA file.
        length_tag: Tag prefix for sequence length in ``S`` lines.

    Returns:
        A NetworkX ``DiGraph`` representing the assembly graph, where
        each node is a ``(seq_id, orientation)`` tuple and each edge
        carries ``match_bases`` and ``distance`` attributes.

    Raises:
        ValueError: If an ``S`` line is missing the ``length_tag``, or
            an ``L`` line references a segment not defined in any ``S`` line.
    """
    G: nx.DiGraph = nx.DiGraph()
    gfa_seq_lengths: dict[str, int] = {}
    l_lines: list[tuple[str, str, str, str, int]] = []

    with gopen(gfa_file, 'rt') as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('S'):
                cols: list[str] = line.split('\t')
                seq_id: str = cols[1]
                ln_match: re.Match[str] | None = re.search(rf'{length_tag}(\d+)', line)
                if ln_match is None:
                    _logger.error(
                        'Cannot find %s tag in S line for sequence %s in GFA file %s',
                        length_tag, repr(seq_id), repr(gfa_file),
                    )
                    raise ValueError(
                        f'Cannot find {length_tag!r} tag in S line '
                        f'for sequence {seq_id!r} in GFA file {gfa_file!r}.'
                    )
                length: int = int(ln_match.group(1))
                gfa_seq_lengths[seq_id] = length
                forward_node: tuple[str, str] = (seq_id, '+')
                reverse_node: tuple[str, str] = (seq_id, '-')
                G.add_node(forward_node, length=length)
                G.add_node(reverse_node, length=length)
            elif line.startswith('L'):
                cols: list[str] = line.split('\t')
                u_seq_id: str = cols[1]
                u_ori: str = cols[2]
                v_seq_id: str = cols[3]
                v_ori: str = cols[4]
                match_bases: int = (
                    int(cols[5][:-1]) if cols[5].endswith('M') else 0
                )
                l_lines.append((u_seq_id, u_ori, v_seq_id, v_ori, match_bases))

    for u_seq_id, u_ori, v_seq_id, v_ori, match_bases in l_lines:
        for seq_id in (u_seq_id, v_seq_id):
            if seq_id not in gfa_seq_lengths:
                _logger.error(
                    'L line references segment %s not defined in any S line in GFA file %s. '
                    'The GFA file appears to be incomplete or malformed.',
                    repr(seq_id), repr(gfa_file),
                )
                raise ValueError(
                    f'L line references segment {seq_id!r} not defined '
                    f'in any S line in GFA file {gfa_file!r}. '
                    f'The GFA file appears to be incomplete or malformed.'
                )
        u: tuple[str, str] = (u_seq_id, u_ori)
        v: tuple[str, str] = (v_seq_id, v_ori)
        if match_bases > G.nodes[u]['length']:
            _logger.warning(
                'match_bases (%d) > node %s length (%d), clamping distance to 0',
                match_bases, repr(u), G.nodes[u]['length'],
            )
        G.add_edge(
            u, v,
            match_bases=match_bases,
            distance=max(0, G.nodes[u]['length'] - match_bases),
        )

    _complete_missing_complementary_edges(G, gfa_file)

    _logger.info(
        'Parsed GFA %r: %d nodes, %d edges',
        gfa_file, G.number_of_nodes(), G.number_of_edges(),
    )
    return G


def _flip_ori(ori: str) -> str:
    """Flip a strand orientation: ``'+'`` becomes ``'-'`` and vice versa.

    Args:
        ori: Strand orientation, must be ``'+'`` or ``'-'``.

    Returns:
        The flipped orientation.
    """
    return '-' if ori == '+' else '+'


def _complete_missing_complementary_edges(G: nx.DiGraph, gfa_file: str) -> None:
    """Add missing complementary edges to a GFA assembly graph.

    In a valid GFA, every ``L`` line should have a complementary link:
    if ``utgA+ -> utgB-`` exists, then ``utgB+ -> utgA-`` must also exist.
    This function detects missing complementary edges, adds them using the
    same ``match_bases`` value, and logs a warning.

    Args:
        G: The assembly ``DiGraph`` to complete.
        gfa_file: Path to the source GFA file (for log messages).
    """
    edges_to_add: list[tuple[tuple[str, str], tuple[str, str], int]] = []

    for u, v, data in G.edges(data=True):
        v_flip: tuple[str, str] = (v[0], _flip_ori(v[1]))
        u_flip: tuple[str, str] = (u[0], _flip_ori(u[1]))
        if not G.has_edge(v_flip, u_flip):
            edges_to_add.append((v_flip, u_flip, data['match_bases']))

    for v_flip, u_flip, match_bases in edges_to_add:
        _logger.warning(
            'Missing complementary edge %s -> %s in GFA %r, adding with match_bases=%d',
            repr(v_flip), repr(u_flip), gfa_file, match_bases,
        )
        G.add_edge(
            v_flip, u_flip,
            match_bases=match_bases,
            distance=max(0, G.nodes[v_flip]['length'] - match_bases),
        )


def parse_yz_agp(
    agp_file: str,
    max_scaffold: int = 113,
) -> list[list[list[str | float]]]:
    """Parse YZ AGP file and extract utg information for the first N scaffolds.

    Processes contig lines (column 5 == ``'W'``) from the AGP file and
    resolves each contig into one or more unitigs (utgs).  Two contig
    naming formats are supported:

    * Single-utg: ``utg_start_end`` or ``utg_start_end_copyN`` — the
      orientation is taken from column 9.
    * Multi-utg: ``utg1_orient_start_end__utg2_orient_start_end`` (double
      underscore separator ``__``) — each utg carries its own
      orientation, and column 9 is always ``'+'``.

    The ``_copy`` suffix (with optional trailing digits) is stripped from
    utg IDs before matching against GFA segments.

    Args:
        agp_file: Path to the YZ AGP file.
        max_scaffold: Maximum number of scaffolds to process.

    Returns:
        A list of scaffolds, where each scaffold is a list of utg
        entries.  Each utg entry is
        ``[utg_id, orientation, center_in_scaffold, center_in_utg]``,
        where ``center_in_scaffold`` and ``center_in_utg`` are floats
        representing the midpoint coordinates (1-based).

    Raises:
        ValueError: If a contig line cannot be parsed.
    """
    scaffolds: list[list[list[str | float]]] = []
    current_scaffold_id: str | None = None
    current_scaffold_utgs: list[list[str | float]] = []
    scaffold_count: int = 0

    with gopen(agp_file, 'rt') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue

            cols: list[str] = line.strip().split('\t')
            if len(cols) < 9:
                continue

            scaffold_id: str = cols[0]
            scaff_start: int = int(cols[1])
            scaff_end: int = int(cols[2])
            comp_type: str = cols[4]

            if scaffold_id != current_scaffold_id:
                if current_scaffold_id is not None and current_scaffold_utgs:
                    scaffolds.append(current_scaffold_utgs)
                    scaffold_count += 1
                    if scaffold_count >= max_scaffold:
                        break
                current_scaffold_id = scaffold_id
                current_scaffold_utgs = []

            if comp_type != 'W':
                continue

            contig_id: str = cols[5]
            contig_ori: str = cols[8]

            if '__' in contig_id:
                parts = contig_id.split('__')
                part_lengths: list[int] = []
                parsed_parts: list[tuple[str, str, int, int]] = []
                for part in parts:
                    clean_part = re.sub(r'_copy\d*\Z', '', part)
                    m = re.match(r'(.+)_([+-])_(\d+)_(\d+)\Z', clean_part)
                    if m is None:
                        raise ValueError(
                            f'Cannot parse multi-utg part {clean_part!r} '
                            f'in contig {contig_id!r}'
                        )
                    utg_id = m.group(1)
                    utg_ori = m.group(2)
                    utg_start = int(m.group(3))
                    utg_end = int(m.group(4))
                    parsed_parts.append((utg_id, utg_ori, utg_start, utg_end))
                    part_lengths.append(utg_end - utg_start + 1)

                total_length = scaff_end - scaff_start + 1
                if sum(part_lengths) != total_length:
                    _logger.warning(
                        'Contig %s: sum of utg lengths (%d) != scaffold span (%d)',
                        repr(contig_id), sum(part_lengths), total_length,
                    )

                cum_start = scaff_start
                for (utg_id, utg_ori, utg_start, utg_end), plen in zip(
                    parsed_parts, part_lengths
                ):
                    utg_scaff_end = cum_start + plen - 1
                    center_scaffold = (cum_start + utg_scaff_end) / 2.0
                    center_utg = (utg_start + utg_end) / 2.0
                    current_scaffold_utgs.append([
                        utg_id, utg_ori, center_scaffold, center_utg,
                    ])
                    cum_start += plen
            else:
                clean_id = re.sub(r'_copy\d*\Z', '', contig_id)
                m = re.match(r'(.+)_(\d+)_(\d+)\Z', clean_id)
                if m is None:
                    raise ValueError(
                        f'Cannot parse single-utg contig {clean_id!r}'
                    )
                utg_id = m.group(1)
                utg_start = int(m.group(2))
                utg_end = int(m.group(3))
                center_scaffold = (scaff_start + scaff_end) / 2.0
                center_utg = (utg_start + utg_end) / 2.0
                current_scaffold_utgs.append([
                    utg_id, contig_ori, center_scaffold, center_utg,
                ])

        if current_scaffold_utgs and scaffold_count < max_scaffold:
            scaffolds.append(current_scaffold_utgs)

    _logger.info(
        'Parsed YZ AGP %r: %d scaffolds, %d total utgs',
        agp_file,
        len(scaffolds),
        sum(len(s) for s in scaffolds),
    )

    for i, scaffold in enumerate(scaffolds):
        _logger.info('Scaffold %d: %d utgs', i + 1, len(scaffold))

    with open('scaffolds_YZ.txt', 'w') as f:
        for i, scaffold in enumerate(scaffolds):
            f.write(f'# Scaffold {i + 1} ({len(scaffold)} utgs)\n')
            for utg in scaffold:
                f.write('\t'.join(str(x) for x in utg) + '\n')
            f.write('\n')
    _logger.info('Debug output written to scaffolds_YZ.txt')

    return scaffolds


def parse_poj_agp(
    agp_file: str,
) -> list[list[list[str | float]]]:
    """Parse POJ AGP file and extract utg information for Chr scaffolds.

    Processes contig lines (column 5 == ``'W'``) from the AGP file and
    resolves each contig into one unitig (utg).  Only scaffolds whose ID
    starts with ``'Chr'`` are processed.  Three contig naming formats are
    supported:

    * Plain: ``utgxxxxxxx`` — columns 7-8 give utg coordinates, orientation
      from column 9.
    * Suffixed: ``utgxxxxxxx_dx`` — the ``_d<N>`` suffix is stripped,
      same as plain format.
    * Ranged: ``utgxxxxxxx:start-end`` — start-end gives utg
      coordinates (1-based), orientation from column 9.

    Args:
        agp_file: Path to the POJ AGP file.

    Returns:
        A list of scaffolds, where each scaffold is a list of utg
        entries.  Each utg entry is
        ``[utg_id, orientation, center_in_scaffold, center_in_utg]``,
        where ``center_in_scaffold`` and ``center_in_utg`` are floats
        representing the midpoint coordinates (1-based).

    Raises:
        ValueError: If a contig line cannot be parsed.
    """
    scaffolds: list[list[list[str | float]]] = []
    current_scaffold_id: str | None = None
    current_scaffold_utgs: list[list[str | float]] = []

    with gopen(agp_file, 'rt') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue

            cols: list[str] = line.strip().split('\t')
            if len(cols) < 9:
                continue

            scaffold_id: str = cols[0]
            if not scaffold_id.startswith('Chr'):
                if current_scaffold_id is not None:
                    if current_scaffold_utgs:
                        scaffolds.append(current_scaffold_utgs)
                    current_scaffold_id = None
                    current_scaffold_utgs = []
                continue

            scaff_start: int = int(cols[1])
            scaff_end: int = int(cols[2])
            comp_type: str = cols[4]

            if scaffold_id != current_scaffold_id:
                if current_scaffold_id is not None and current_scaffold_utgs:
                    scaffolds.append(current_scaffold_utgs)
                current_scaffold_id = scaffold_id
                current_scaffold_utgs = []

            if comp_type != 'W':
                continue

            contig_id: str = cols[5]
            contig_start: int = int(cols[6])
            contig_end: int = int(cols[7])
            contig_ori: str = cols[8]

            if ':' in contig_id:
                m = re.match(r'(.+):(\d+)-(\d+)\Z', contig_id)
                if m is None:
                    raise ValueError(
                        f'Cannot parse ranged contig {contig_id!r}'
                    )
                utg_id = m.group(1)
                utg_start = int(m.group(2))
                utg_end = int(m.group(3))
            else:
                clean_id = re.sub(r'_d\d+\Z', '', contig_id)
                utg_id = clean_id
                utg_start = contig_start
                utg_end = contig_end

            center_scaffold = (scaff_start + scaff_end) / 2.0
            center_utg = (utg_start + utg_end) / 2.0
            current_scaffold_utgs.append([
                utg_id, contig_ori, center_scaffold, center_utg,
            ])

        if current_scaffold_utgs:
            scaffolds.append(current_scaffold_utgs)

    _logger.info(
        'Parsed POJ AGP %r: %d scaffolds, %d total utgs',
        agp_file,
        len(scaffolds),
        sum(len(s) for s in scaffolds),
    )

    for i, scaffold in enumerate(scaffolds):
        _logger.info('Scaffold %d: %d utgs', i + 1, len(scaffold))

    with open('scaffolds_POJ.txt', 'w') as f:
        for i, scaffold in enumerate(scaffolds):
            f.write(f'# Scaffold {i + 1} ({len(scaffold)} utgs)\n')
            for utg in scaffold:
                f.write('\t'.join(str(x) for x in utg) + '\n')
            f.write('\n')
    _logger.info('Debug output written to scaffolds_POJ.txt')

    return scaffolds


def compare_graph_scaffolds(
    scaffolds: list[list[list[str | float]]],
    graph: nx.DiGraph,
    label: str,
) -> dict[str, int]:
    """Compare scaffolds adjacency against the assembly graph.

    For each pair of adjacent utgs in each scaffold, checks whether
    a path exists in ``graph`` for each of the four orientation
    combinations (``++``, ``+-``, ``-+``, ``--``), and classifies
    the result into one of several categories.

    Categories:

    * **not_connected**: none of the four combinations has a path.
    * **wrong_order**: the scaffold's combination has no path, but at
      least one other combination does.
    * **direct_edge**: the scaffold's combination is a direct edge in
      ``graph``.
    * **path_le_250kb**: the scaffold's combination has a path whose
      total ``distance`` ≤ 250 Kb, but is not a direct edge.
    * **path_250kb_500kb**: total ``distance`` ∈ (250 Kb, 500 Kb].
    * **path_500kb_750kb**: total ``distance`` ∈ (500 Kb, 750 Kb].
    * **path_750kb_1mb**: total ``distance`` ∈ (750 Kb, 1 Mb].

    Args:
        scaffolds: Scaffolds parsed from AGP, each scaffold is a list
            of ``[utg_id, orientation, center_in_scaffold, center_in_utg]``.
        graph: Assembly graph built from GFA, with ``distance`` edge
            attribute for shortest-path computation.
        label: Label used in log messages (e.g. ``'YZ'`` or ``'POJ'``).

    Returns:
        A dict with keys: ``total``, ``not_connected``,
        ``wrong_order``, ``direct_edge``, ``path_le_250kb``,
        ``path_250kb_500kb``, ``path_500kb_750kb``, ``path_750kb_1mb``,
        ``center_diffs`` (list of ``(float, str, str, str, str, str)``
        tuples, each tuple contains the scaffold-minus-graph center
        distance, the connection category (``'direct_edge'``,
        ``'has_path'``, or ``'wrong_order'``), and the two utg IDs
        and orientations (``id_a``, ``ori_a``, ``id_b``, ``ori_b``),
        for all pairs that have at least one working orientation
        combo in ``graph``).
    """
    total = 0
    not_connected = 0
    wrong_order = 0
    direct_edge = 0
    path_le_250kb = 0
    path_250kb_500kb = 0
    path_500kb_750kb = 0
    path_750kb_1mb = 0
    center_diffs: list[tuple[float, str, str, str, str, str]] = []

    for scaffold in scaffolds:
        for i in range(len(scaffold) - 1):
            total += 1
            id_a: str = str(scaffold[i][0])
            ori_a: str = str(scaffold[i][1])
            id_b: str = str(scaffold[i + 1][0])
            ori_b: str = str(scaffold[i + 1][1])

            scaffold_combo = ((id_a, ori_a), (id_b, ori_b))

            all_combos: list[tuple[tuple[str, str], tuple[str, str]]] = []
            for oa in ('+', '-'):
                for ob in ('+', '-'):
                    all_combos.append(((id_a, oa), (id_b, ob)))

            path_exists: dict[
                tuple[tuple[str, str], tuple[str, str]], bool
            ] = {}
            for src, tgt in all_combos:
                if src in graph and tgt in graph:
                    path_exists[(src, tgt)] = nx.has_path(graph, src, tgt)
                else:
                    path_exists[(src, tgt)] = False

            scaffold_has_path = path_exists[scaffold_combo]
            other_combos = [c for c in all_combos if c != scaffold_combo]
            any_other_has_path = any(path_exists[c] for c in other_combos)

            if not scaffold_has_path and not any_other_has_path:
                not_connected += 1
                continue

            center_scaff_a: float = float(scaffold[i][2])
            center_scaff_b: float = float(scaffold[i + 1][2])
            center_utg_a: float = float(scaffold[i][3])
            center_utg_b: float = float(scaffold[i + 1][3])
            scaffold_center_dist: float = abs(
                center_scaff_b - center_scaff_a
            )

            graph_combo: tuple[tuple[str, str], tuple[str, str]]
            raw_path_len: float

            category: str
            if not scaffold_has_path and any_other_has_path:
                best_raw_len: float = float('inf')
                graph_combo = other_combos[0]
                for c in other_combos:
                    if not path_exists[c]:
                        continue
                    c_src, c_tgt = c
                    raw_len = (
                        nx.shortest_path_length(
                            graph, c_src, c_tgt, weight='distance',
                        ) - graph.nodes[c_src]['length']
                    )
                    if raw_len < best_raw_len:
                        best_raw_len = raw_len
                        graph_combo = c
                raw_path_len = best_raw_len
                if best_raw_len > 1_000_000:
                    not_connected += 1
                    continue
                wrong_order += 1
                category = 'wrong_order'
            else:
                graph_combo = scaffold_combo
                src, tgt = graph_combo

                if graph.has_edge(src, tgt):
                    direct_edge += 1
                    category = 'direct_edge'
                    raw_path_len = (
                        graph.edges[src, tgt]['distance']
                        - graph.nodes[src]['length']
                    )
                else:
                    raw_path_len = (
                        nx.shortest_path_length(
                            graph, src, tgt, weight='distance',
                        ) - graph.nodes[src]['length']
                    )
                    if raw_path_len > 1_000_000:
                        not_connected += 1
                        continue
                    category = 'has_path'
                    path_len: float = max(0.0, raw_path_len)
                    if path_len <= 250_000:
                        path_le_250kb += 1
                    elif path_len <= 500_000:
                        path_250kb_500kb += 1
                    elif path_len <= 750_000:
                        path_500kb_750kb += 1
                    else:
                        path_750kb_1mb += 1

            src, tgt = graph_combo
            ori_graph_a: str = graph_combo[0][1]
            ori_graph_b: str = graph_combo[1][1]

            if ori_graph_a == '+':
                exit_dist: float = (
                    graph.nodes[src]['length'] - center_utg_a
                )
            else:
                exit_dist: float = center_utg_a - 1

            if ori_graph_b == '+':
                entry_dist: float = center_utg_b
            else:
                entry_dist: float = (
                    graph.nodes[tgt]['length'] - center_utg_b + 1
                )

            graph_center_dist: float = (
                exit_dist + raw_path_len + entry_dist
            )
            center_diffs.append(
                (scaffold_center_dist - graph_center_dist, category,
                 id_a, ori_a, id_b, ori_b)
            )

    if center_diffs:
        diff_values: list[float] = [d[0] for d in center_diffs]
        mean_diff: float = statistics.mean(diff_values)
        median_diff: float = statistics.median(diff_values)
        stdev_diff: float = (
            statistics.stdev(diff_values)
            if len(diff_values) > 1
            else 0.0
        )
        _logger.info(
            '=== %s Center-Distance Diff (scaffold - graph) ===',
            label,
        )
        _logger.info(
            '  n=%d mean=%.1f median=%.1f stdev=%.1f',
            len(diff_values), mean_diff, median_diff, stdev_diff,
        )
        _logger.info(
            '  min=%.1f max=%.1f',
            min(diff_values), max(diff_values),
        )

    denom = total if total else 1
    _logger.info('=== %s Graph-Scaffold Comparison ===', label)
    _logger.info('  Total adjacent utg pairs: %d', total)
    _logger.info(
        '  Not connected in graph: %d (%.2f%%)',
        not_connected, 100.0 * not_connected / denom,
    )
    _logger.info(
        '  Wrong order/orientation: %d (%.2f%%)',
        wrong_order, 100.0 * wrong_order / denom,
    )
    _logger.info(
        '  Direct edge: %d (%.2f%%)',
        direct_edge, 100.0 * direct_edge / denom,
    )
    _logger.info(
        '  Path <= 250 Kb: %d (%.2f%%)',
        path_le_250kb, 100.0 * path_le_250kb / denom,
    )
    _logger.info(
        '  Path (250 Kb, 500 Kb]: %d (%.2f%%)',
        path_250kb_500kb, 100.0 * path_250kb_500kb / denom,
    )
    _logger.info(
        '  Path (500 Kb, 750 Kb]: %d (%.2f%%)',
        path_500kb_750kb, 100.0 * path_500kb_750kb / denom,
    )
    _logger.info(
        '  Path (750 Kb, 1 Mb]: %d (%.2f%%)',
        path_750kb_1mb, 100.0 * path_750kb_1mb / denom,
    )

    return {
        'total': total,
        'not_connected': not_connected,
        'wrong_order': wrong_order,
        'direct_edge': direct_edge,
        'path_le_250kb': path_le_250kb,
        'path_250kb_500kb': path_250kb_500kb,
        'path_500kb_750kb': path_500kb_750kb,
        'path_750kb_1mb': path_750kb_1mb,
        'center_diffs': center_diffs,
    }


def write_category_counts(
    yz_result: dict,
    poj_result: dict,
    output_path: str,
) -> None:
    """Write connection category counts to a CSV file for pie chart plotting in R.

    Args:
        yz_result: Result dict from ``compare_graph_scaffolds`` for YZ.
        poj_result: Result dict from ``compare_graph_scaffolds`` for POJ.
        output_path: Path to the output CSV file.
    """
    level1_categories: list[str] = [
        'not_connected',
        'wrong_order',
        'direct_edge',
    ]
    level2_categories: list[str] = [
        'path_le_250kb',
        'path_250kb_500kb',
        'path_500kb_750kb',
        'path_750kb_1mb',
    ]
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('genome,category,level,count,percentage\n')
        for genome, result in [('YZ', yz_result), ('POJ', poj_result)]:
            total = result['total']
            denom = total if total else 1
            for cat in level1_categories:
                cnt = result[cat]
                pct = round(100.0 * cnt / denom, 2)
                f.write(f'{genome},{cat},1,{cnt},{pct}\n')
            has_path_cnt = sum(result[c] for c in level2_categories)
            has_path_pct = round(100.0 * has_path_cnt / denom, 2)
            f.write(f'{genome},has_path,1,{has_path_cnt},{has_path_pct}\n')
            has_path_denom = has_path_cnt if has_path_cnt else 1
            for cat in level2_categories:
                cnt = result[cat]
                pct = round(100.0 * cnt / has_path_denom, 2)
                f.write(f'{genome},{cat},2,{cnt},{pct}\n')
    _logger.info('Category counts written to %s', repr(output_path))


def write_center_diffs(
    yz_result: dict,
    poj_result: dict,
    output_path: str,
) -> None:
    """Write center_diffs data to a CSV file for distribution plotting in R.

    Args:
        yz_result: Result dict from ``compare_graph_scaffolds`` for YZ.
        poj_result: Result dict from ``compare_graph_scaffolds`` for POJ.
        output_path: Path to the output CSV file.
    """
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('genome,center_diff,category,id_a,ori_a,id_b,ori_b\n')
        for genome, result in [('YZ', yz_result), ('POJ', poj_result)]:
            for diff_val, cat, id_a, ori_a, id_b, ori_b in result['center_diffs']:
                f.write(f'{genome},{diff_val},{cat},{id_a},{ori_a},{id_b},{ori_b}\n')
    _logger.info('Center diffs written to %s', repr(output_path))


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
    )
    args: Namespace = parse_args()
    check_file_exists(args.yz_gfa, args.yz_agp, args.poj_gfa, args.poj_agp)
    validate_input_files(args)

    yz_graph: nx.DiGraph = parse_gfa(args.yz_gfa)
    poj_graph: nx.DiGraph = parse_gfa(args.poj_gfa)
    yz_scaffolds = parse_yz_agp(args.yz_agp)
    poj_scaffolds = parse_poj_agp(args.poj_agp)

    yz_result = compare_graph_scaffolds(yz_scaffolds, yz_graph, 'YZ')
    poj_result = compare_graph_scaffolds(poj_scaffolds, poj_graph, 'POJ')

    write_category_counts(yz_result, poj_result, 'category_counts.csv')
    write_center_diffs(yz_result, poj_result, 'center_diffs.csv')


if __name__ == '__main__':
    main()

