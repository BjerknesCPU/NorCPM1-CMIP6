#!/usr/bin/env perl
#-----------------------------------------------------------------------------------------------
#
# lt_archive.pl
#
# This utility moves files to the long term archive 
#
#-----------------------------------------------------------------------------------------------
use strict;
#use warnings;
#use diagnostics;

use Cwd;
use English;
use Getopt::Long;
use IO::File;
use IO::Handle;
use File::Find;
#-----------------------------------------------------------------------------------------------

sub usage {
    die <<EOF;
SYNOPSIS
     lt_archive [options]
OPTIONS
     User supplied values are denoted in angle brackets (<>).  Any value that contains
     white-space must be quoted.  Long option names may be supplied with either single
     or double leading dashes.  A consequence of this is that single letter options may
     NOT be bundled.

     -help [or -h]        Print usage to STDOUT.
     -silent [or -s]      Turns on silent mode - only fatal messages issued.
     -verbose [or -v]     Turn on verbose echoing of settings made by configure.
EOF
}

#-----------------------------------------------------------------------------------------------
# Setting autoflush (an IO::Handle method) on STDOUT helps in debugging.  It forces the test
# descriptions to be printed to STDOUT before the error messages start.

*STDOUT->autoflush();                  

#-----------------------------------------------------------------------------------------------
# Set the directory that contains the CCSM configuration scripts.  If the create_newcase command was
# issued using a relative or absolute path, that path is in $ProgDir.  Otherwise assume the
# command was issued from the current working directory.

(my $ProgName = $0) =~ s!(.*)/!!;      # name of this script
my $ProgDir = $1;                      # name of directory containing this script -- may be a
                                       # relative or absolute path, or null if the script is in
                                       # the user's PATH
my $cwd = getcwd();                    # current working directory
my $cfgdir;                            # absolute pathname of directory that contains this script
if ($ProgDir) { 
    $cfgdir = absolute_path($ProgDir);
    chdir("$cfgdir/..");
    $cwd = getcwd();
    print "cwd = $cwd\n";
} else {
    $cfgdir = $cwd;
}

#-----------------------------------------------------------------------------------------------
# Parse command-line options.
my %opts = (
	    cache       => "config_cache.xml",
	    );
GetOptions(
    "h|help"                    => \$opts{'help'},
    "s|silent"                  => \$opts{'silent'},
    "v|verbose"                 => \$opts{'verbose'},
)  or usage();

# Give usage message.
usage() if $opts{'help'};

# Check for unparsed argumentss
if (@ARGV) {
    print "ERROR: unrecognized arguments: @ARGV\n";
    usage();
}

(-f "env_case.xml")  or  die <<"EOF";
** Cannot find env_case.xml **
EOF

# Define 3 print levels:
# 0 - only issue fatal error messages
# 1 - only informs what files are created (default)
# 2 - verbose
my $print = 1;
if ($opts{'silent'})  { $print = 0; }
if ($opts{'verbose'}) { $print = 2; }
my $eol = "\n";

my %cfg = ();           # build configuration

#-----------------------------------------------------------------------------------------------
# Make sure we can find required perl modules and configuration files.
# Look for them in the directory that contains the configure script.

# Check for the configuration definition file.
my $config_def_file = "config_definition.xml";
(-f "$cfgdir/$config_def_file")  or  die <<"EOF";
** Cannot find configuration definition file \"$config_def_file\" in directory ./Tools **
EOF

# The XML::Lite module is required to parse the XML configuration files.
(-f "$cfgdir/XML/Lite.pm")  or  die <<"EOF";
** Cannot find perl module \"XML/Lite.pm\" in directory ./Tools **
EOF

# The ConfigCase module provides utilities to store and manipulate the configuration.
(-f "$cfgdir/ConfigCase.pm")  or  die <<"EOF";
** Cannot find perl module \"ConfigCase.pm\" in directory ./Tools **
EOF

if ($print>=2) { print "Setting configuration directory to $cfgdir$eol"; }

#-----------------------------------------------------------------------------------------------
# Add $cfgdir/perl5lib to the list of paths that Perl searches for modules
my @dirs = (  $cfgdir, "$cfgdir/Tools");
unshift @INC, @dirs;
require XML::Lite;
require ConfigCase;

my $config_file = "$cfgdir/$config_def_file";

my $cfg_ref = ConfigCase->new("$config_file", "env_case.xml");
$cfg_ref->reset_setup("env_build.xml");
$cfg_ref->reset_setup("env_conf.xml");
$cfg_ref->reset_setup("env_run.xml");
$cfg_ref->reset_setup("env_mach_pes.xml");


my $dout_s_root = $cfg_ref->getresolved('DOUT_S_ROOT');
my $dout_l_msroot = $cfg_ref->getresolved('DOUT_L_MSROOT');


# this bit of code handles the bluefire users who
# suffered the transition from mss to hpss

my (@pfiles,@dfiles);


my $tstamp = `date +%Y%m%d%H%M%S%N`;
chomp $tstamp;  # removes trailing return character

my $newdir = $dout_s_root."/.lta$tstamp";

my $accnt = $cfg_ref->get('DOUT_L_HPSS_ACCNT');

my $cnt = hpss_put($dout_s_root, $newdir, $accnt);

rmtree($newdir) if($cnt<=0);

opendir(S,".");
my @oldltadirs = grep(/.lta/,readdir(S));
closedir(S);

foreach $newdir (@oldltadirs){
    my $age = (-A $newdir);
    if($age > 0.5){       # age in days of the directory
	my $cnt = hpss_put($dout_s_root, $newdir, $accnt);
	rmtree($newdir) if($cnt<=0);
    }
}

#-----------------------------------------------------------------------------------------------
# FINISHED ####################################################################################
#-----------------------------------------------------------------------------------------------



sub plainfiles{
    my ($dev,$ino,$mode,$nlink,$uid,$gid);
    (($dev,$ino,$mode,$nlink,$uid,$gid) = lstat($_)) &&
	-f _ && push(@pfiles, $File::Find::name);
}

sub dirfiles{
    my ($dev,$ino,$mode,$nlink,$uid,$gid);
    (($dev,$ino,$mode,$nlink,$uid,$gid) = lstat($_)) &&
	-d _ && push(@dfiles, $File::Find::name);
    
}

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}

sub rmtree{
    my ($dir) = @_;

    finddepth(\&dirfiles, 'find', $dir);
    foreach(@dfiles){
	print "Removing $_\n" if($print>1);
	rmdir($_) || warn "rmtree $_: $!\n";
    }
    $#dfiles=-1;
}

sub hpss_put{
    my ($olddir, $newdir, $accnt) = @_;
    my $cnt=0;

    mkdir $newdir;
    print "mv $olddir/* $newdir" if($print>1);
    system("mv $olddir/* $newdir");

    find(\&plainfiles, $newdir);
    $cnt = $#pfiles;
    $#pfiles=-1;
    chdir $newdir;

    my $hsicmd = "hsi ";
    $hsicmd.="-a $accnt" if($accnt > 0);
    $hsicmd .=" \"mkdir -p -m 0775 $dout_l_msroot; cd $dout_l_msroot; put -dPR *; chmod -R +r *\"";
    
    print "$hsicmd\n" if($print>0);
    system($hsicmd) if($cnt>0);

    chdir $dout_s_root;

    find(\&plainfiles, $newdir);
    $cnt = $#pfiles;
    print "Files still in $newdir: @pfiles \n" if($#pfiles>=0 && $print>0);
    $#pfiles=-1;
}
    
sub absolute_path {
#
# Convert a pathname into an absolute pathname, expanding any . or .. characters.
# Assumes pathnames refer to a local filesystem.
# Assumes the directory separator is "/".
#
  my $path = shift;
  my $cwd = getcwd();  # current working directory
  my $abspath;         # resulting absolute pathname

# Strip off any leading or trailing whitespace.  (This pattern won't match if
# there's embedded whitespace.
  $path =~ s!^\s*(\S*)\s*$!$1!;

# Convert relative to absolute path.

  if ($path =~ m!^\.$!) {          # path is "."
      return $cwd;
  } elsif ($path =~ m!^\./!) {     # path starts with "./"
      $path =~ s!^\.!$cwd!;
  } elsif ($path =~ m!^\.\.$!) {   # path is ".."
      $path = "$cwd/..";
  } elsif ($path =~ m!^\.\./!) {   # path starts with "../"
      $path = "$cwd/$path";
  } elsif ($path =~ m!^[^/]!) {    # path starts with non-slash character
      $path = "$cwd/$path";
  }

  my ($dir, @dirs2);
  my @dirs = split "/", $path, -1;   # The -1 prevents split from stripping trailing nulls
                                     # This enables correct processing of the input "/".

  # Remove any "" that are not leading.
  for (my $i=0; $i<=$#dirs; ++$i) {
      if ($i == 0 or $dirs[$i] ne "") {
	  push @dirs2, $dirs[$i];
      }
  }
  @dirs = ();

  # Remove any "."
  foreach $dir (@dirs2) {
      unless ($dir eq ".") {
	  push @dirs, $dir;
      }
  }
  @dirs2 = ();

  # Remove the "subdir/.." parts.
  foreach $dir (@dirs) {
    if ( $dir !~ /^\.\.$/ ) {
        push @dirs2, $dir;
    } else {
        pop @dirs2;   # remove previous dir when current dir is ..
    }
  }
  if ($#dirs2 == 0 and $dirs2[0] eq "") { return "/"; }
  $abspath = join '/', @dirs2;
  return( $abspath );
}


