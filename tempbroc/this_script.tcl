# VMD TCL Script to Display All Bond Lengths and Angles
# Usage: vmd -dispdev text -e this_script.tcl > output.txt
#    or: source this_script.tcl (in VMD console)

# Load your PDB file
set pdbfile "dfhbi_CG.pdb"
mol new $pdbfile

# Get the molecule ID
set molid [molinfo top]

# Manually define bonds based on your connectivity
# Bonds: A7A2 A2A3 A2A1 A1A3 A3A4 A4A6 A4A5 A5A6
# Note: VMD uses 0-based indexing, so atom 1 (B1) = index 0

puts "\nManually defining bonds based on specified connectivity..."

# Define the bond list
# Format: {atom1_index atom2_index}
set manual_bonds {
    {6 1}  
    {1 2}  
    {1 0}  
    {0 2}  
    {2 3}  
    {3 5}  
    {3 4}  
    {4 5}  
}

# Build bonding information for each atom
set natoms 7
array set atom_bonds {}
for {set i 0} {$i < $natoms} {incr i} {
    set atom_bonds($i) {}
}

# Populate the bond arrays
foreach bond $manual_bonds {
    set atom1 [lindex $bond 0]
    set atom2 [lindex $bond 1]
    lappend atom_bonds($atom1) $atom2
    lappend atom_bonds($atom2) $atom1
}

# Apply bonds to each atom
for {set i 0} {$i < $natoms} {incr i} {
    set sel [atomselect top "index $i"]
    $sel setbonds [list $atom_bonds($i)]
    $sel delete
}

puts "Bonds defined: B7-B2, B2-B3, B2-B1, B1-B3, B3-B4, B4-B6, B4-B5, B5-B6"

# Function to calculate distance between two atoms
proc calc_distance {atom1 atom2} {
    set coord1 [lindex [[atomselect top "index $atom1"] get {x y z}] 0]
    set coord2 [lindex [[atomselect top "index $atom2"] get {x y z}] 0]
    return [veclength [vecsub $coord2 $coord1]]
}

# Function to calculate angle between three atoms
proc calc_angle {atom1 atom2 atom3} {
    set coord1 [lindex [[atomselect top "index $atom1"] get {x y z}] 0]
    set coord2 [lindex [[atomselect top "index $atom2"] get {x y z}] 0]
    set coord3 [lindex [[atomselect top "index $atom3"] get {x y z}] 0]
    
    set vec1 [vecsub $coord1 $coord2]
    set vec2 [vecsub $coord3 $coord2]
    
    set dot [vecdot $vec1 $vec2]
    set len1 [veclength $vec1]
    set len2 [veclength $vec2]
    
    if {$len1 == 0 || $len2 == 0} {
        return 0
    }
    
    set cosang [expr $dot / ($len1 * $len2)]
    # Clamp to [-1, 1] to avoid numerical errors
    if {$cosang > 1.0} {set cosang 1.0}
    if {$cosang < -1.0} {set cosang -1.0}
    
    return [expr {acos($cosang) * 180.0 / 3.14159265359}]
}

# Get all atoms
set all [atomselect top all]
set natoms [$all num]

puts "\n========================================="
puts "STRUCTURE INFORMATION"
puts "=========================================\n"
puts "Total atoms: $natoms"
puts "Residue: [lindex [$all get resname] 0]"
puts "Residue ID: [lindex [$all get resid] 0]"

puts "\n========================================="
puts "BOND LENGTH ANALYSIS"
puts "=========================================\n"

# Analyze all bonds
set bondlist {}
set bondcount 0
for {set i 0} {$i < $natoms} {incr i} {
    set sel [atomselect top "index $i"]
    set bonds [lindex [$sel getbonds] 0]
    set name1 [lindex [$sel get name] 0]
    set resid1 [lindex [$sel get resid] 0]
    set resname1 [lindex [$sel get resname] 0]
    set index1 [lindex [$sel get index] 0]
    
    foreach j $bonds {
        if {$j > $i} {  # Only count each bond once
            set sel2 [atomselect top "index $j"]
            set name2 [lindex [$sel2 get name] 0]
            set resid2 [lindex [$sel2 get resid] 0]
            set resname2 [lindex [$sel2 get resname] 0]
            set index2 [lindex [$sel2 get index] 0]
            
            set dist [calc_distance $i $j]
            
            incr bondcount
            puts [format "Bond %3d: %4s (atom %2d) -- %4s (atom %2d)  Length: %.4f Å" \
                $bondcount $name1 [expr $index1+1] $name2 [expr $index2+1] $dist]
            
            lappend bondlist [list $i $j]
            $sel2 delete
        }
    }
    $sel delete
}

if {$bondcount == 0} {
    puts "\nWARNING: No bonds detected!"
    puts "VMD's automatic bond detection may not have found bonds."
    puts "This could mean:"
    puts "  1. Atoms are too far apart for standard bonding"
    puts "  2. You need to manually define bonds"
    puts "  3. You need a PSF/topology file"
}

puts "\n========================================="
puts "ANGLE ANALYSIS"
puts "=========================================\n"

# Analyze all angles
set anglecount 0
for {set i 0} {$i < $natoms} {incr i} {
    set sel [atomselect top "index $i"]
    set bonds [lindex [$sel getbonds] 0]
    
    # For each central atom with at least 2 bonds
    if {[llength $bonds] >= 2} {
        set name2 [lindex [$sel get name] 0]
        set resid2 [lindex [$sel get resid] 0]
        set resname2 [lindex [$sel get resname] 0]
        set index2 [lindex [$sel get index] 0]
        
        # Get all pairs of bonded atoms
        for {set m 0} {$m < [llength $bonds]} {incr m} {
            for {set n [expr $m + 1]} {$n < [llength $bonds]} {incr n} {
                set atom1 [lindex $bonds $m]
                set atom3 [lindex $bonds $n]
                
                set sel1 [atomselect top "index $atom1"]
                set sel3 [atomselect top "index $atom3"]
                
                set name1 [lindex [$sel1 get name] 0]
                set resid1 [lindex [$sel1 get resid] 0]
                set resname1 [lindex [$sel1 get resname] 0]
                set index1 [lindex [$sel1 get index] 0]
                
                set name3 [lindex [$sel3 get name] 0]
                set resid3 [lindex [$sel3 get resid] 0]
                set resname3 [lindex [$sel3 get resname] 0]
                set index3 [lindex [$sel3 get index] 0]
                
                set angle [calc_angle $atom1 $i $atom3]
                
                incr anglecount
                puts [format "Angle %3d: %4s (atom %2d) -- %4s (atom %2d) -- %4s (atom %2d)  Angle: %.2f°" \
                    $anglecount \
                    $name1 [expr $index1+1] \
                    $name2 [expr $index2+1] \
                    $name3 [expr $index3+1] \
                    $angle]
                
                $sel1 delete
                $sel3 delete
            }
        }
    }
    $sel delete
}

puts "\n========================================="
puts "SUMMARY"
puts "=========================================\n"
puts "Total atoms: $natoms"
puts "Total bonds detected: $bondcount"
puts "Total angles calculated: $anglecount"

# Print atomic coordinates for reference
puts "\n========================================="
puts "ATOMIC COORDINATES"
puts "=========================================\n"
set all [atomselect top all]
for {set i 0} {$i < $natoms} {incr i} {
    set sel [atomselect top "index $i"]
    set name [lindex [$sel get name] 0]
    set coord [lindex [$sel get {x y z}] 0]
    set x [lindex $coord 0]
    set y [lindex $coord 1]
    set z [lindex $coord 2]
    puts [format "Atom %2d (%4s): x=%8.3f y=%8.3f z=%8.3f" [expr $i+1] $name $x $y $z]
    $sel delete
}

puts "\n"
$all delete

# Optional: quit VMD after analysis (comment out if running in GUI)
# quit