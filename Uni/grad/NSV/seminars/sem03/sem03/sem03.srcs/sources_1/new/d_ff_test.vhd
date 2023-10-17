library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity d_ff_test is
    Port ( clk : in STD_LOGIC;
           d : in STD_LOGIC;
           q : out STD_LOGIC);
end d_ff_test;

architecture Behavioral of d_ff_test is

--signal t: std_logic := '0';

begin

process(clk)
--defining t as a variable instead of a signal results in a completely different schematic
variable t: std_logic := '0';
begin
    if rising_edge(clk) then
        t := d;
        q <= t;
    end if;
end process;

end Behavioral;
