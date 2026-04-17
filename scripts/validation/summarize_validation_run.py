#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path


LOG_STAGE_NAMES = {
    "01_filter.log": "filter-vcf",
    "02_prepare.log": "prepare",
    "03_compute.log": "compute",
    "04_plot_static.log": "plot_static",
    "05_plot_interactive.log": "plot_interactive",
}


def human_size(num_bytes):
    value = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024.0 or unit == "TB":
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{num_bytes} B"


def parse_log(log_path):
    text = log_path.read_text(encoding="utf-8", errors="replace")
    info = {
        "stage": LOG_STAGE_NAMES.get(log_path.name, log_path.stem),
        "path": str(log_path),
        "completed": False,
        "real_seconds": None,
        "user_seconds": None,
        "sys_seconds": None,
        "kept_records": None,
        "total_records": None,
        "text": text,
    }

    real_match = re.search(r"^real\s+([0-9.]+)$", text, re.MULTILINE)
    user_match = re.search(r"^user\s+([0-9.]+)$", text, re.MULTILINE)
    sys_match = re.search(r"^sys\s+([0-9.]+)$", text, re.MULTILINE)
    count_match = re.search(
        r"\((\d+) kept / (\d+) total records\)",
        text,
    )

    if real_match:
        info["completed"] = True
        info["real_seconds"] = float(real_match.group(1))
    if user_match:
        info["user_seconds"] = float(user_match.group(1))
    if sys_match:
        info["sys_seconds"] = float(sys_match.group(1))
    if count_match:
        info["kept_records"] = int(count_match.group(1))
        info["total_records"] = int(count_match.group(2))

    return info


def collect_artifacts(workflow_root):
    artifacts = []
    for path in sorted(workflow_root.iterdir()):
        if path.name == "logs":
            continue
        if path.name == ".DS_Store":
            continue
        if path.is_file():
            artifacts.append(
                {
                    "name": path.name,
                    "path": str(path),
                    "type": "file",
                    "size_bytes": path.stat().st_size,
                    "size_human": human_size(path.stat().st_size),
                }
            )
        elif path.is_dir():
            size_bytes = sum(
                child.stat().st_size
                for child in path.rglob("*")
                if child.is_file()
            )
            artifacts.append(
                {
                    "name": path.name,
                    "path": str(path),
                    "type": "directory",
                    "size_bytes": int(size_bytes),
                    "size_human": human_size(size_bytes),
                }
            )
    return artifacts


def build_markdown(args, stage_logs, artifacts):
    lines = [
        f"# {args.title}",
        "",
        f"Workflow root: `{args.workflow_root}`",
        "",
        f"Provenance mode: `{args.provenance_mode}`",
        "",
        "## Summary",
        "",
        args.summary,
        "",
        "## Stages",
        "",
    ]

    for info in stage_logs:
        lines.append(f"### {info['stage']}")
        lines.append("")
        lines.append(f"- log: `{info['path']}`")
        lines.append(f"- completed: `{info['completed']}`")
        if info["real_seconds"] is not None:
            lines.append(f"- wall time: `{info['real_seconds']:.2f} s`")
        if info["user_seconds"] is not None:
            lines.append(f"- user CPU time: `{info['user_seconds']:.2f} s`")
        if info["sys_seconds"] is not None:
            lines.append(f"- system CPU time: `{info['sys_seconds']:.2f} s`")
        if info["kept_records"] is not None and info["total_records"] is not None:
            lines.append(
                f"- retained records: `{info['kept_records']:,}` of `{info['total_records']:,}`"
            )
        lines.append("")

    lines.extend(
        [
            "## Artifacts",
            "",
        ]
    )

    for artifact in artifacts:
        lines.append(
            f"- `{artifact['name']}`: `{artifact['type']}`, `{artifact['size_human']}`"
        )

    if args.machine_specs:
        lines.extend(
            [
                "",
                "## Machine",
                "",
                f"Machine specs source: `{args.machine_specs}`",
            ]
        )

    if args.script_path:
        lines.extend(
            [
                "",
                "## Source script",
                "",
                f"`{args.script_path}`",
            ]
        )

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Summarize a local validation workflow into JSON and Markdown."
    )
    parser.add_argument("workflow_root", type=Path)
    parser.add_argument("--title", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--provenance-mode", required=True)
    parser.add_argument("--script-path")
    parser.add_argument("--machine-specs")
    parser.add_argument("--out-prefix", default="validation_summary")
    args = parser.parse_args()

    workflow_root = args.workflow_root.resolve()
    logs_dir = workflow_root / "logs"
    log_paths = [logs_dir / name for name in LOG_STAGE_NAMES if (logs_dir / name).exists()]
    stage_logs = [parse_log(path) for path in log_paths]
    artifacts = collect_artifacts(workflow_root)

    summary = {
        "title": args.title,
        "workflow_root": str(workflow_root),
        "provenance_mode": args.provenance_mode,
        "summary": args.summary,
        "script_path": args.script_path,
        "machine_specs": args.machine_specs,
        "stages": [
            {k: v for k, v in stage.items() if k != "text"}
            for stage in stage_logs
        ],
        "artifacts": artifacts,
    }

    out_json = workflow_root / f"{args.out_prefix}.json"
    out_md = workflow_root / f"{args.out_prefix}.md"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    out_md.write_text(build_markdown(args, stage_logs, artifacts), encoding="utf-8")


if __name__ == "__main__":
    main()
