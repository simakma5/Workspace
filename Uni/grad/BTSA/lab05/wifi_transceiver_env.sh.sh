#!/bin/bash
# Bash Menu Script Example

if [ "$EUID" -ne 0 ]
then echo "Please run this script as root"
  exit 1
fi


PS3="Please select your classroom table row number 1-5:"
options=("Row 1. - (close to the whiteboard)" "Row 2." "Row 3." "Row 4." "Row 5.(close to the door)" "Quit.")

select opt in "${options[@]}"
do

    case $opt in

	    "Row 1. - (close to the whiteboard)" | "Row 2." | "Row 3."  | "Row 4." | "Row 5.(close to the door)" )
            ROW=$REPLY
	    echo "You have selected: $opt"
	    break
	    ;;
	"Quit.")
	    echo "Terminating!"
	    exit 1
            break
            ;;
        *)
                echo "invalid option $REPLY !"
                ;;
    esac
done
PS3="Please select table position:"
options=("Left (window side)" "Right (door side)" "Quit.")

select opt in "${options[@]}"
do

    case $opt in

	 "Left (window side)")
        SIDE=2
	    OTHERSIDE=1
	    echo "You have selected: $opt"
	    break
	    ;;
    "Right (door side)")
	    SIDE=1
	    OTHERSIDE=2
	    echo"You have selected: $opt"
	    break
	    ;;
	"Quit.")
	    echo "Terminating!"
	    exit 1
            break
            ;;
        *)
                echo "invalid option $REPLY !"
                ;;
    esac
done

echo "...deleting possible previously made TAP device"
ip tuntap del tap0 mode tap || { echo "$LINENO line command failed" ; }

echo "...adding tuntap device tap0 in mode: tap"
ip tuntap add dev tap0 mode tap || { echo "$LINENO line command failed" ; exit 1; }

echo "...turning the tap0 device DOWN"
ip link set dev tap0 down || { echo "$LINENO line command failed" ; exit 1; }


echo "...setting the mac address of our tap device to:  12:34:56:78:90:$ROW$SIDE"
ip link set dev tap0 address 12:34:56:78:90:$ROW$SIDE || { echo "$LINENO line command failed" ; exit 1; } # replace ab by a = row number b = 1 for left PC, b = 2 for right PC

echo "...setting the MTU"
ip link set dev tap0 mtu 440 || { echo "$LINENO line command failed" ; exit 1; }

echo "...setting the device to UP"
ip link set dev tap0 up || { echo "$LINENO line command failed" ; exit 1; }

echo "...setting the IP address of tap0 to :  192.168.$ROW.$SIDE/24"
ip addr add 192.168.$ROW.$SIDE/24 dev tap0 || { echo "$LINENO line command failed" ; exit 1; } # ROW_NUMBER = 120 + row number, PC_IP is 1 for left PC and 2 for right PC


echo "...deleting possible previous route"
ip route del  192.168.$ROW.0/24 || { echo "$LINENO line command failed, but not critical" ; } # Note that route may not exists previously

echo "...adding the new route with MSS setting: "
ip route add 192.168.$ROW.0/24 dev tap0 advmss 400 || { echo "$LINENO line command failed" ; exit 1; }

#echo "... deleting possible previous tc qdisc netem delay if any was present"
#tc qdisc del dev tap0 root || { echo "$LINENO line command failed, but not critical" ; } # Note that device may not exists previously


#replace may work better than add here
echo "...adding tc qdisc netem delay 10ms"
tc qdisc replace dev tap0 root netem delay 10ms || { echo "$LINENO line command failed" ; exit 1; }


echo "adding the ARP record: arp -s 192.168.$ROW.$OTHERSIDE 12:34:56:78:90:$ROW$OTHERSIDE"
arp -s 192.168.$ROW.$OTHERSIDE 12:34:56:78:90:$ROW$OTHERSIDE || { echo "$LINENO line command failed" ; exit 1; } # replace ab by a = row number !! b = 1 for right PC, b = 2 for left PC !!

exit 0
