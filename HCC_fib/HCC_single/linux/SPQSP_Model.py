import os
import sys
import subprocess
import xml.etree.ElementTree as ET
from shutil import copyfile
from datetime import datetime

class SPQSP_Model:
    def __init__(self, execFile, paramFileName='../resource/param_all_test.xml', outputDir='/home/szhan121/scratch4-apopel1/szhan121/MARCC/HCC/pyABC_output', time=280, grid=1, grid_interval=1):
        self.ParamFileName = paramFileName
        self.seed = 1
        self.OutputDir = outputDir
        self.Time = time
        self.executable = os.path.abspath(execFile)
        self.grid = grid
        self.grid_interval = grid_interval
        self.outputs_folder = outputDir  # Folder to store the output files
        if os.name == 'nt':
            self.executable = self.executable.replace(os.altsep, os.sep) + '.exe'  # Windows path handling

    def RunModel(self, SampleID='', Parameters={}):
        # Get the output path
        self.set_outputPath(SampleID)

        # Create a copy of the parameter file and modify it if parameters are provided
        xml_file_out = os.path.join(self.OutputDir, f'param_modified_{SampleID}.xml')
        original_param_file = self.ParamFileName  # Use original parameter file

        if len(Parameters) > 0:
            self.generate_xml_file(original_param_file, xml_file_out, Parameters)
        else:
            # If no parameters are provided, just copy the original file
            copyfile(original_param_file, xml_file_out)

        # Update the parameter file for the simulation
        self.ParamFileName = xml_file_out
        print(self.ParamFileName)
        
        # Execute the simulation with the modified parameter file
        callingModel = [
            self.executable, '-s', str(self.seed),
            '-p', self.ParamFileName, '-o', self.OutputDir,
            '-t', str(self.Time), '-G', str(self.grid),
            '--grid-interval', str(self.grid_interval)
        ]
        
        cache = subprocess.run(callingModel, universal_newlines=True, capture_output=True)
        if cache.returncode != 0:
            print(f"Error: model output error! Executable: {self.executable}. returned: \n{str(cache.returncode)}")
            print(f"Error message (stderr): {cache.stderr}")
            return -1
        else:
            print(cache.stdout)
            return 0
        
    def set_outputPath(self, sampleID):
        output_path_base = '/home/szhan121/scratch4-apopel1/szhan121/MARCC/HCC/pyABC_output'
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.OutputDir = os.path.join(output_path_base, f"output_{sampleID}_{timestamp}")
        #self.OutputDir = f"{self.OutputDir}_{sampleID}"
        os.makedirs(self.OutputDir, exist_ok=True)  # Ensure the full folder is created
        
    def get_outputPath(self):
        return self.OutputDir

    def generate_xml_file(self, xml_file_in, xml_file_out, dic_parameters):
        # Copy the original XML file to the new location
        copyfile(xml_file_in, xml_file_out)
        tree = ET.parse(xml_file_out)
        xml_root = tree.getroot()
        print("xml_root: ", xml_root)
        # Modify parameters in the copied XML file based on the dictionary keys
        for key, val in dic_parameters.items():
            # Split the key into parent and child (e.g., "MDSC.moveProb" becomes "MDSC" and "moveProb")
            if '.' in key:
                parent_element, child_element = key.split('.')
                set_xml_elements_value(xml_root, child_element, val, parent_element=parent_element)
            else:
                # For other parameters that aren't nested, just update them directly
                set_xml_elements_value(xml_root, key, val)

        # Save the modified XML tree to the output file
        tree.write(xml_file_out)
    

def set_xml_elements_value(xml_root, element_name, new_value, parent_element=None):
    """
    Find and set the value of a given XML element by name, optionally under a specific parent element.

    Parameters:
    - root: The root element of the parsed XML tree.
    - element_name: The name of the XML element to find.
    - new_value: The new value to set for the element.
    - parent_element: The name of the parent element to restrict the search to (e.g., 'MDSC').
    """
    xml_mapping = {
        "MDSC_moveProb": "MDSC.moveProb",
        "k_rec_MDSC": "k_rec_MDSC",
        "CancerCell_moveSteps": "CancerCell.moveSteps",
        "TCD4_moveSteps": "TCD4.moveSteps",
        "TCell_moveSteps": "TCell.moveSteps",
        "Fib_moveSteps": "Fib.moveSteps",
        "Mac_moveSteps": "Mac.moveSteps"
        # Add other mappings as needed
    }

    # Use the mapping to convert element_name back to the XML format
    if element_name in xml_mapping:
        element_name = xml_mapping[element_name]

    if parent_element:
        # Find the parent element first (e.g., 'MDSC')
        parent = xml_root.find(f".//{parent_element}")
        if parent is not None:
            # Now find the child element (e.g., 'moveProb') within the parent element
            element = parent.find(f".//{element_name.replace('.', '/')}")
        else:
            print(f"Parent element '{parent_element}' not found.")
            return
    else:
        # Find the element without specifying a parent (global search)
        element = xml_root.find(f".//{element_name.replace('.', '/')}")
    
    if element is not None:
        print(f"Original value of {parent_element}.{element_name}: {element.text}")
        element.text = str(new_value[0])  # Update the value
        print(f"New value of {parent_element}.{element_name}: {element.text}")
    else:
        print(f"Element '{element_name}' not found under '{parent_element}' in the XML.")

def get_xml_element_value(xml_root, key):
    elem = xml_root.findall(key)
    if (len(elem) != 1):
        raise ValueError(f"""
                Multiples occurrences or none found to this key: {key}, occurences: {[pos.text for pos in elem]}.
                Key examples:
                # key cell cycle example: ".//*[@name='CD8 Tcell']/phenotype/cycle/phase_transition_rates/rate[4]"
                # key substrates example: ".//*[@name='TNF']/physical_parameter_set/diffusion_coefficient"
                # key parameter example: ".//random_seed"
                """)
    return elem[0].text

def check_parameter_in_xml(xml_file_in, key_parameter):
    tree = ET.parse(xml_file_in)
    xml_root = tree.getroot()
    try: text_elem = get_xml_element_value(xml_root, key_parameter)
    except ValueError as e:
        print(f"Error in Parameters definition: {e}")
        sys.exit(1)
