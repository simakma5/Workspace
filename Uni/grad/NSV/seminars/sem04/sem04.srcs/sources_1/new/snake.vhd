library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity snake is
    Port ( clk : in STD_LOGIC;
           rst_n : in STD_LOGIC;
           dir : in STD_LOGIC;
           en : in STD_LOGIC;
           led : out STD_LOGIC_VECTOR (15 downto 0));
end snake;

architecture Behavioral of snake is
    -- Shift register initialization
    signal sreg_i: std_logic_vector(15 downto 0) := X"0001";

    begin
        led <= sreg_i;
        -- trigger process on every clk change
        process(clk) begin
            -- preform the following code on clk rising edge
            if rising_edge(clk) then
                if rst_n = '0' then
                    sreg_i <= X"0001";
                else
                    if dir = '1' and en = '1' then
                        sreg_i <= sreg_i(14 downto 0) & sreg_i(15);
                    elsif dir = '0' and en = '1' then
                        sreg_i <= sreg_i(0) & sreg_i(15 downto 1);
                    else
                        sreg_i <= sreg_i;
                    end if;
                end if;
            end if;
        end process;


end Behavioral;
