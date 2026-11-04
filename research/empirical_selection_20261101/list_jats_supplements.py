import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("path", type=Path)
args = parser.parse_args()
tree = ET.parse(args.path)
for node in tree.iter():
    if node.tag.rsplit("}", 1)[-1] == "supplementary-material":
        links = [value for child in node.iter() for key, value in child.attrib.items() if key.rsplit("}", 1)[-1] == "href"]
        print(json.dumps({"id": node.attrib.get("id"), "text": " ".join(node.itertext()).strip(), "links": links}))
